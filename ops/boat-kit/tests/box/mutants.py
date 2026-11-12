"""KIT10 mutants: one exact replacement (it must occur once) in a copy of the candidate, then the test that should
fail on it. KIT10 round 2 (REVIEW-KIT10 L5): each mutant names the RESULT id it must fail, and it counts as killed
only when the test exits non-zero and that id's line says DIFFERS; a mutant whose test fails elsewhere is
WRONG-ID, not killed. JOBS at a time, at most half the vCPUs; the install mutants run after the rest, one at a time,
since install-kit10.sh refuses while any bash runs drive.sh, which the drive mutants' tests do. Writes one line per
mutant and a summary to OUT, and exits 1 unless every mutant is killed. KIT10_MUT_KINDS (a comma list of kinds, as
runner,h1box) runs only those, for a check on the controller; the chain sets none. KIT10 round 3 (REVIEW-KIT10-R2 L1)
adds three install mutants, one of them on INSTALL.md's rollback (kind installmd).
usage: python3 -I -B mutants.py HARNESS WORK OUT [JOBS]"""
import concurrent.futures as cf
import os
import shutil
import subprocess
import sys
import time

H, WK, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
HALF = max(1, (os.cpu_count() or 2) // 2)
JOBS = min(int(sys.argv[4]), HALF) if len(sys.argv) > 4 else HALF
RUNNER, DRIVE, KIT, T = H + "/cand/grm-boat.sh", H + "/cand/drive.sh", H + "/kit", H + "/tests"
DRAFTS, HBK = H + "/drafts", H + "/handback"
MARK_OLD = '''      if [ -d "$A/out/data/species/$s" ]; then
        cp -al "$A/out/data/species/$s" "$SD/out/data/species/" || die 2 "could not carry $s's outputs from $A"
      fi
      ln "$f" "$SD/out/run/gates/" || die 2 "could not carry $s's done mark from $A"
'''
MARK_NEW = '''      ln "$f" "$SD/out/run/gates/" || die 2 "could not carry $s's done mark from $A"
      if [ -d "$A/out/data/species/$s" ]; then
        cp -al "$A/out/data/species/$s" "$SD/out/data/species/" || die 2 "could not carry $s's outputs from $A"
      fi
'''
M = [
    # runner (test-runner-kit10.sh)
    ("R-stage-shared", "runner", '  mkdir -p "$base" && STAGE=$(mktemp -d "$base/call-$1-XXXXXXXX") || fail', '  mkdir -p "$base" && STAGE=$base || fail', "S1-parallel-push-data-same-lane"),
    ("R-unstage-removes-lane", "runner", 'unstage() { rm -rf "$STAGE"; return 0; }', 'unstage() { rm -rf "$STAGE"; rmdir "${STAGE%/*}" 2> /dev/null; return 0; }', "S3b-lane-folder-kept"),
    ("R-window-from-launch", "runner", "o pid rc t0 tw s s0=0 need started=0 blind=0 gap", "o pid rc t0 tw s s0=0 need started=1 blind=0 gap", "G1-slow-start-then-healthy-completes"),
    ("R-blind-ignored", "runner", '[ $blind -ge "$SCP_BLIND" ] && SCP_WHY=', 'false && SCP_WHY=', "U1-unreadable-size-stops"),
    ("R-no-start-bound", "runner", '      elif [ $(( $(date +%s) - t0 )) -ge "$SCP_START" ]; then', '      elif false; then', "P4-stopped-well-before-the-timeout"),
    ("R-floor-0-still-on", "runner", '    [ "$SCP_ON" = 1 ] || continue\n', '', "P5-timeout-rc-124"),
    ("R-extend-may-shorten", "runner", '  [ "$want" -gt "$a0" ] || fail 2', '  true || fail 2', "X1-shorter-refused-unchanged"),
    ("R-extend-no-auto-stop", "runner", 'out=$(bt extend "$ID" --ttl "$ttl" 2>&1)', 'out=$(bt extend "$ID" --ttl "$ttl" --no-auto-stop 2>&1)', "X2-longer-set"),
    ("R-extend-no-after-check", "runner", '  [ "$a1" -ge "$a0" ] || fail', '  true || fail', "X9-moved-back-within-120s-86"),
    ("R-extend-short-ok", "runner", '  [ "$a1" -ge $(( want - 120 )) ] || fail', '  true || fail', "X6-extend-that-changed-nothing-86"),
    ("R-extend-window-wide", "runner", '  [ "$a1" -ge $(( want - 120 )) ] || fail', '  [ "$a1" -ge $(( want - 900 )) ] || fail', "X7-extend-600s-short-86"),
    ("R-extend-window-zero", "runner", '  [ "$a1" -ge $(( want - 120 )) ] || fail', '  [ "$a1" -ge "$want" ] || fail', "X8-extend-60s-short-within-120s-0"),
    ("R-extend-no-end-line", "runner", 'echo "END $a1"', ':', "X2b-end-line-is-the-new-end"),
    ("R-vocab-any-bytes", "runner", "VOCAB.get(parts[1]) != h", "parts[1] not in VOCAB", "V2-real-pins-refuse-other-bytes"),
    ("R-vocab-partial-ok", "runner", "if vocab and vocab != sorted(VOCAB):", "if False:", "V4-partial-refused"),
    ("R-vocab-name-anywhere", "runner", '    if parts[-1].lower() in VOCAB:', '    if False:', "W4-tables-sigdrug-name-refused"),
    ("R-vocab-name-case", "runner", 'parts[-1].lower() in VOCAB', 'parts[-1] in VOCAB', "W5-sigdrug-name-in-other-case-refused"),
    ("R-combine-any-dest", "runner", 'combine_dest = dest == "/home/user/grm-combine"', 'combine_dest = True', "W8-same-root-to-another-destination-refused"),
    ("R-allowed-ignores-dest", "runner", 'or (combine_dest and combine.fullmatch(p)))\ndef is_receipt', 'or combine.fullmatch(p))\ndef is_receipt', "W9-master-manifest-at-the-top-elsewhere-refused"),
    ("R-h1work-outside-parts", "runner", r'r"parts/%s/h1work/(%s/){0,5}%s\.(json|jsonl|tsv|txt|md)"', r'r"(%s/){0,3}h1work/(%s/){0,5}%s\.(json|jsonl|tsv|txt|md)"', "W10-h1work-outside-parts-refused-on-combine"),
    ("R-h1work-any-file", "runner", r'r"parts/%s/h1work/(%s/){0,5}%s\.(json|jsonl|tsv|txt|md)"', r'r"parts/%s/h1work/(%s/){0,5}%s\.[a-z]+"', "H2b-text-member-with-another-suffix-refused"),
    ("R-master-any-depth", "runner", r'               r"MASTER-MANIFEST[A-Za-z0-9._-]*\.sha256"]', r'               r"(?:[A-Za-z0-9][A-Za-z0-9._-]*/){0,3}MASTER-MANIFEST[A-Za-z0-9._-]*\.sha256"]', "W11-master-manifest-in-a-part-refused-on-combine"),
    ("R-no-mem-retry", "runner", '    sleep "${GRM_MEM_RETRY_SECS:-5}"; win_mem; read -r free commit <<< "$M"', '    sleep "${GRM_MEM_RETRY_SECS:-5}"; read -r free commit <<< "$M"', "M-empty-then-ok"),
    # kit chain-lib.sh (test-runner-kit10.sh's C checks)
    ("C-finish-passes-a-failure", "chain", '"$CHAIN_LOG/STEPS"; exit 1', '"$CHAIN_LOG/STEPS"; exit 0', "C1-failed-middle-step-fails-the-chain"),
    ("C-must-goes-on", "chain", 'must() { step "$@" || { echo "$(date -u +%FT%TZ) chain stopped after $1" >> "$CHAIN_LOG/STEPS"; finish; }; }', 'must() { step "$@"; }', "C3-must-stops-the-chain"),
    ("C-finish-no-step-passes", "chain", '  if [ "$CHAIN_N" = 0 ]; then', '  if false; then', "C5-no-step-fails"),
    # kit h1-box.sh and r11-box.sh (test-h1-box.sh)
    ("B-union-count-unchecked", "h1box", '[ "$n" -gt 0 ] && [ "$n" = "$u" ] ||', 'true ||', "U2-union-count-differs-stops"),
    ("B-fasta-set-not-union-ok", "h1box", '= "$(cat "$R/union.txt")" ] ||', '= "$(cat "$R/union.txt")" ] || true ||', "V2-fasta-set-not-the-union"),
    ("B-fasta-set-extra-ok", "h1box", "grep -v '^fasta/')\" ] ||", "grep -v '^fasta/')\" ] || true ||", "V3-fasta-set-holds-more"),
    ("B-raw-writable", "h1box", '  chmod -R a-w "$RAW"\n', '  :\n', "V1-verify"),
    ("B-estimate-not-checked", "h1box", 'sys.exit(0 if d.get("ok") is True and d.get("union_without_fasta") == 0 else 1)', 'sys.exit(0)', "V4-estimate-not-ok"),
    ("B-gate-line-not-needed", "h1box", "grep -q 'revision-12 gate: PASSED' \"$R/estimate.err\" ||", "true ||", "V5-real-spec-without-the-gate-line"),
    ("B-marker-any-spec", "h1box", '"$SPECF" ||\n    RC=2 fail "the metadata holds receipts/', '"$SPECF" || true ||\n    RC=2 fail "the metadata holds receipts/', "M1-marker-on-a-real-spec-refused"),
    ("B-extra-any", "h1box", '[ "$H1_EXTRA" = --below-panel-floor ] ||', 'true ||', "E1b-other-extra-on-a-synthetic-spec-refused"),
    ("B-run-on-writable-raw", "h1box", '[ -d "$RAW" ] && [ ! -w "$RAW" ] ||', 'true ||', "R3-writable-raw-refused"),
    ("B-vocab-set-extra-ok", "h1box", '= vocabulary ] ||', '= vocabulary ] || true ||', "U4-vocabulary-set-holds-more"),
    ("B-set-not-rehashed", "h1box", '(cd "$d" && sha256sum -c --strict --quiet .grm-input-manifest.sha256) ||', 'true ||', "U5-pushed-file-differs-from-its-manifest"),
    ("B-estimate-at-raw-free", "h1box", '--estimate-only --free-disk-bytes "$use"', '--estimate-only --free-disk-bytes "$free"', "V7-estimate-at-usable-free"),
    ("B-run-at-raw-free", "h1box", '--jobs "${H1_JOBS:-2}" --free-disk-bytes "$use"', '--jobs "${H1_JOBS:-2}" --free-disk-bytes "$free"', "R4-run-at-usable-free"),
    ("B-usable-no-floor", "h1box", '[[ "$u" =~ ^[0-9]+$ ]] && [ "$u" -gt 0 ] ||', 'true ||', "V8-free-below-the-pause-floor-4"),
    ("B-marker-move-unchecked", "h1box", 'rmdir "$1/receipts" || RC=2 fail', 'rmdir "$1/receipts" || true || RC=2 fail', "M3-marker-that-cannot-move-2"),
    ("K-status-never-paused", "r11", '      if [ -n "$gp" ]; then', '      if false; then', "D4-guard-pause-reads-paused"),
    # kit h1-manifest.py and the drive (test-drive-h1.sh)
    ("F-missing-fasta-ok", "h1man", "        if miss:\n            die(", "        if False:\n            die(", "H8-union-genome-without-fasta-2-box-closed-reason-named"),
    ("F-meta-takes-features", "h1man", 'or p.split("/")[0] in ("catalog", "species-check")]', 'or p.split("/")[0] in ("catalog", "species-check", "bvbrc")]', "H1-meta-set-no-features-curated"),
    ("D-no-reference-wait", "drive", '&& die 79 "the H-1 reference box', '&& false && die 79 "the H-1 reference box', "H5-run-box-before-the-reference-79"),
    ("D-vocab-pin-unchecked", "drive", '[ "$VSHA" = "$VOCAB_MANIFEST_SHA256" ] ||', 'true ||', "H7-vocabulary-pin-mismatch-2-no-box"),
    ("D-needs-one-file", "drive", 'for n in ${NEEDS//,/ }; do [ "$n" != - ] && [ ! -e "$n" ] && die 79 "needs $n"; done', '[ "$NEEDS" != - ] && [ ! -e "$NEEDS" ] && die 79 "needs $NEEDS"', "H4-combine-box-rc"),
    ("D-ttl-column-ignored", "drive", '[ "$KIND" = h1 ] && TTL=$(awk -v h="$TTL_H"', 'false && TTL=$(awk -v h="$TTL_H"', "H1-ttl-from-column-and-template"),
    ("D-floor-not-exported", "drive", '  [ -n "${GRM_SCP_FLOOR_MBPS:-}" ] && export GRM_SCP_FLOOR_MBPS', '  :', "H1-floor-5-reaches-the-runner"),
    ("D-extend-raw-boat", "drive", 'GRM_BOAT_ADOPT=1 "$R" extend "$ID" --ttl "$rem" > "$SD/extend.out" 2>&1', 'bt extend "$ID" --ttl "$rem" > "$SD/extend.out" 2>&1', "H14-extend-through-the-runner"),
    ("D-src-commit-unchecked", "drive", '[ "$C" = "$SRC_COMMIT" ] || die 2 "$SRC is at $C, not SRC_COMMIT', 'true || die 2 "$SRC is at $C, not SRC_COMMIT', "H9-other-source-commit-2-no-box"),
    ("D-combine-from-template", "drive", '[ "$H1ROLE" = combine ] && NEWTMPL=""', ':', "H4-combine-default-no-template"),
    ('D-combine-flags-unchecked', 'drive', '[ "$(cut -d\' \' -f2- <<< "$pf" | sort -u | wc -l)" = 1 ] && ! grep -q \': none$\' <<< "$pf" ||', 'true ||', "H4b-combine-parts-with-other-flags-2-no-box"),
    ("D-flags-shape-unchecked", "drive", '''[ -f "$f" ] && [ "$(wc -l < "$f")" = 1 ] && grep -qxE 'H1_TRACE=[01] H1_EXTRA=[^ ]*' "$f"''', '[ -f "$f" ]', "H4d-combine-parts-with-empty-flags-2-no-box"),
    ("D-deadline-asked-not-got", "drive", '    DEADLINE=$end;', '    DEADLINE=$(( T0 + want ));', "H14b-deadline-is-the-end-the-box-got"),
    ("D-extend-end-optional", "drive", '[ -n "$end" ] || bail 86', '[ -n "$end" ] || end=$(( T0 + want )); true || bail 86', "H14c-extend-without-end-86-closed"),
    ("D-paused-waits", "drive", '  if [ "$state" = paused ]; then', '  if false; then', "H17-guard-pause-stops-81-with-its-reason"),
    ("D-git-status-rc-ignored", "drive", '[ $src_rc = 0 ] || die 2 "git status failed', 'true || die 2 "git status failed', "H18-git-status-failing-2-no-box"),
    ("D-row-placeholder-ok", "drive", '[ -z "$ph" ] || die 2 "stage $BOX: placeholder left in $STAGES', 'true || die 2 "stage $BOX: placeholder left in $STAGES', "H19-stage-table-placeholder-2-no-box"),
    # the drive's REVIEW-KIT9 L1 and L2, KIT10's L3 and the H-1 dry run's stale still-up entry (test-drive-l2-kit10.sh)
    ("D-copy2-no-delete", "drivel2", 'rsync -a --delete "$SD/out/" "$COPY2/$BOX/"', 'rsync -a "$SD/out/" "$COPY2/$BOX/"', "L1-stale-gone-from-copy2"),
    ("D-mark-first", "drivel2", MARK_OLD, MARK_NEW, "L2-no-done-mark-without-outputs"),
    ('D-resume-dirs-left-read-only', 'drivel2', '  [ -d "$RES" ] && find "$RES" -type d -exec chmod u+rwx {} \\;\n', '', "L3-read-only-resume-set-rc"),
    ('D-close-keeps-stale-alive', 'drivel2', '    x="$ALIVE "; while [[ "$x" == *" $ID "* ]]; do x=${x/ "$ID" / }; done; ALIVE=${x% }\n', '', "L4-gone-box-not-86"),
    # install-kit10.sh, REVIEW-KIT10 L4 (test-install-kit10.sh)
    ("I-no-dir-write-check", "install", '  [ -d "${dst%/*}" ] && [ -w "${dst%/*}" ] ||', '  true ||', "I7-read-only-target-dir-stops-before-any-change"),
    ("I-temp-file-left", "install", '  [ -n "$t" ] && rm -f "$t"\n', '', "I9-failure-after-mktemp-leaves-no-temp"),
    ("I-new-files-left", "install", '    if [ "$old" = new ]; then rm -f "$dst" || ok=0; continue; fi\n', '    if [ "$old" = new ]; then continue; fi\n', "I10-last-move-fails-all-seven-undone"),
    ("I-pre-kit10-kept", "install", '    rm -rf "$P"\n', '', "I8-step-2-failure-rolls-back"),
    ("I-pgrep-any-mention", "install", "if pgrep -af '^(/[^ ]*/)?bash( -[A-Za-z]+)* ([^ ]*/)?(drive|grm-boat|chain-tests)\\.sh( |$)'; then",
     "if pgrep -af 'drive\\.sh|grm-boat\\.sh|chain-tests\\.sh'; then", "I12-a-command-naming-drive-sh-does-not-stop-it"),
    # install-kit10.sh and INSTALL.md, REVIEW-KIT10-R2 L1 (test-install-kit10.sh)
    ("I-no-signal-trap", "install", 'trap \'stop2 "interrupted by a signal"\' INT TERM HUP\n', '', "I14-term-in-step-2-rolls-back"),
    ("I-listed-after-the-mv", "install", '  IN+=("$src|$dst|$old|$new|$ref")\n  mv "$t" "$dst" || stop2 "could not move $t to $dst"\n  t=""\n',
     '  mv "$t" "$dst" || stop2 "could not move $t to $dst"\n  t=""; IN+=("$src|$dst|$old|$new|$ref")\n', "I17-term-just-after-a-move-rolls-that-file-back"),
    ("I-md-keeps-temp-files", "installmd", '    rm -f $D/drive/drive.sh.kit10-* $R/grm-boat.sh.kit10-* $D/kit/*.kit10-*\n', '',
     "I18-install-md-rollback-removes-left-temp-files"),
]
KINDS = [k for k in os.environ.get("KIT10_MUT_KINDS", "").split(",") if k]
if KINDS:
    M = [m for m in M if m[1] in KINDS]
FILE = {"runner": "grm-boat.sh", "chain": "kit/chain-lib.sh", "h1box": "kit/h1-box.sh", "r11": "kit/r11-box.sh",
        "h1man": "kit/h1-manifest.py", "drive": "drive.sh", "drivel2": "drive.sh", "install": "handback/install-kit10.sh",
        "installmd": "handback/INSTALL.md"}


def verdict(out, want):
    """KILLED when want's RESULT line says DIFFERS; else WRONG-ID (another check failed) or NO-LINE."""
    seen = False
    for l in out.splitlines():
        f = l.split()
        if len(f) > 1 and f[0] == "RESULT" and f[1] == want:
            seen = True
            if " DIFFERS" in l:
                return "KILLED", l
    return ("WRONG-ID" if seen else "NO-LINE"), ""


def one(m):
    mid, kind, old, new, want = m
    d = os.path.join(WK, mid)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    kit = os.path.join(d, "kit")
    shutil.copytree(KIT, kit)
    runner, drive, hb = os.path.join(d, "grm-boat.sh"), os.path.join(d, "drive.sh"), os.path.join(d, "handback")
    shutil.copy2(RUNNER, runner); shutil.copy2(DRIVE, drive)
    if kind in ("install", "installmd"):
        shutil.copytree(HBK, hb)
    target = os.path.join(d, FILE[kind])
    t = open(target).read()
    if t.count(old) != 1:
        return mid, "NOT-APPLIED", 0, want, "%d matches" % t.count(old)
    open(target, "w").write(t.replace(old, new))
    cmd = {"runner": ["bash", T + "/test-runner-kit10.sh", runner, d + "/w", kit + "/chain-lib.sh"],
           "chain": ["bash", T + "/test-runner-kit10.sh", RUNNER, d + "/w", kit + "/chain-lib.sh"],
           "h1box": ["bash", T + "/test-h1-box.sh", d + "/w", kit],
           "r11": ["bash", T + "/test-h1-box.sh", d + "/w", kit],
           "h1man": ["bash", T + "/test-drive-h1.sh", d + "/w", DRIVE, kit, DRAFTS],
           "drive": ["bash", T + "/test-drive-h1.sh", d + "/w", drive, KIT, DRAFTS],
           "drivel2": ["bash", T + "/test-drive-l2-kit10.sh", d + "/w", drive],
           "install": ["bash", T + "/test-install-kit10.sh", d + "/w", hb, H],
           "installmd": ["bash", T + "/test-install-kit10.sh", d + "/w", hb, H]}[kind]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    open(os.path.join(d, "test.out"), "w").write(p.stdout + p.stderr)
    v, line = verdict(p.stdout, want)
    if p.returncode == 0:
        v = "SURVIVED"
    if v != "KILLED":
        line = next((l for l in p.stdout.splitlines() if " DIFFERS" in l), "")
    return mid, v, round(time.time() - t0), want, line[:150]


os.makedirs(WK, exist_ok=True)
print("mutants %d, jobs %d (vCPUs %s)" % (len(M), JOBS, os.cpu_count()), flush=True)
res = []
with cf.ThreadPoolExecutor(JOBS) as ex:
    for r in ex.map(one, [m for m in M if m[1] not in ("install", "installmd")]):
        res.append(r)
        print("%s %s %ss want %s | %s" % r, flush=True)
for m in M:
    if m[1] in ("install", "installmd"):
        r = one(m)
        res.append(r)
        print("%s %s %ss want %s | %s" % r, flush=True)
k = sum(1 for r in res if r[1] == "KILLED")
with open(OUT, "w") as f:
    f.write("mutants %d, jobs %d (vCPUs %s)\n" % (len(M), JOBS, os.cpu_count()))
    for r in res:
        f.write("%s %s %ss want %s | %s\n" % r)
    f.write("MUTANTS %d killed of %d (killed = the named RESULT id DIFFERS)\n" % (k, len(res)))
print("MUTANTS %d killed of %d (killed = the named RESULT id DIFFERS)" % (k, len(res)))
sys.exit(0 if k == len(res) else 1)
