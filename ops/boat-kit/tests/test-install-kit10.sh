#!/usr/bin/env bash
# KIT10 round 2: install-kit10.sh on a fake laptop tree built from the harness's live copies (no laptop, every file
# synthetic or a harness copy). One "RESULT id want got verdict detail" line per check, then "FAILS n"; exits 1 on
# any failure. I1-I6 are round 1's probe/test-install.sh cases (clean install, INSTALL.md's rollback, install after
# it, second run, drifted live file, changed handback file). I7-I12 are REVIEW-KIT10 L4: its cases A (I1), B (a
# read-only target dir, I7) and C (a step after mktemp fails, I9), a failure after the first file is in (I8), a
# failed move of the last file (I10), and the running-script check (I11, I12). Round 3: I14-I18 are REVIEW-KIT10-R2
# L1, a signal in step 2 and INSTALL.md's rollback with a temp file left.
# usage: test-install-kit10.sh WORK HANDBACK LIVE
#   HANDBACK: a dir with install-kit10.sh, INSTALL.md and next/; LIVE: work/k10 (live/drive.sh, live/grm-boat.sh and
#   PORT's kit under port/kit, the copies the install's live hashes name)
set -u
WK=${1:?work} HBS=${2:?handback} LW=${3:?live work dir}
chmod -R u+w "$WK" 2> /dev/null; rm -rf "$WK"; mkdir -p "$WK"
WK=$(cd "$WK" && pwd)
FAILS=0
res() { local v=as-predicted; [ "$2" = "$3" ] || { v=DIFFERS; FAILS=$((FAILS + 1)); }; printf 'RESULT %s want=%s got=%s %s %s\n' "$1" "$2" "$3" "$v" "${4:-}"; }
yn() { if "$@"; then echo yes; else echo no; fi; }
[ "$(id -u)" != 0 ] || echo "NOTE running as root: a read-only dir doesn't stop root, so I7 can't hold"
O=$WK/F.orig F=$WK/F HB=$WK/HB
mkdir -p "$O/port/drive" "$O/port/kit" "$O/runner"
cp "$LW/live/drive.sh" "$O/port/drive/drive.sh"; cp "$LW/live/grm-boat.sh" "$O/runner/grm-boat.sh"
for f in mtb-box.sh r11-box.sh stage-manifest.py ckpt.sh guard.sh pg-run.sh; do cp "$LW/port/kit/$f" "$O/port/kit/$f"; done
chmod 755 "$O/port/drive/drive.sh" "$O/runner/grm-boat.sh" "$O/port/kit/"*.sh; chmod 644 "$O/port/kit/stage-manifest.py"
git -C "$O/port/kit" init -q; git -C "$O/port/kit" add .; git -C "$O/port/kit" -c user.name=t -c user.email=t@t commit -qm live
mkdir -p "$HB"; cp -a "$HBS/install-kit10.sh" "$HBS/INSTALL.md" "$HBS/next" "$HB/"
(cd "$HB" && find . -type f ! -name SHA256SUMS | LC_ALL=C sort | xargs sha256sum > SHA256SUMS)
reset() { chmod -R u+w "$F" 2> /dev/null; rm -rf "$F" "$HB/pre-kit10"; cp -a "$O" "$F"; }
tree() { # every entry of the fake tree (but the kit's .git) with its mode and, for a file, its sha256
  (cd "$1" && find . -path ./port/kit/.git -prune -o -print | LC_ALL=C sort | while read -r p; do
    if [ -f "$p" ]; then echo "$(stat -c %a "$p") $(sha256sum < "$p" | cut -c1-16) $p"; else echo "$(stat -c %a "$p") - $p"; fi; done); }
same() { [ "$(tree "$O")" = "$(tree "$F")" ]; }
inst() { (cd "${2:-$HB}" && KIT10_D=$F/port KIT10_R=$F/runner PATH=${SHIM:+$SHIM:}$PG:$PATH bash install-kit10.sh) > "$WK/out.$1" 2>&1; }
temps() { find "$F" -name '*.kit10-*' | wc -l; }
push() { git -C "$F/port/kit" ls-files -co --exclude-standard | LC_ALL=C sort | paste -sd' '; }
last() { tail -n 1 "$WK/out.$1" | cut -c1-140; }
# the chain that runs this test is itself "bash chain-tests.sh": a pgrep shim drops this test's own ancestors from
# the real pgrep's lines, so the install's pattern still sees every other process
PG=$WK/pg; mkdir -p "$PG"; REALPG=$(command -v pgrep)
cat > "$PG/pgrep" << EOF
#!/usr/bin/env bash
anc=" "; p=\$PPID; while [ -n "\$p" ] && [ "\$p" -gt 1 ]; do anc="\$anc\$p "; p=\$(ps -o ppid= -p "\$p" | tr -d ' '); done
"$REALPG" "\$@" | while read -r pid rest; do [[ "\$anc" == *" \$pid "* ]] || echo "\$pid \$rest"; done | grep . && exit 0; exit 1
EOF
chmod +x "$PG/pgrep"
# INSTALL.md's rollback block, with the two roots pointed at the fake tree
awk '/^    D=\/home\/layth\//{f=1} f && /^    /{print substr($0, 5); next} f{exit}' "$HB/INSTALL.md" |
  sed "s#/home/layth/review/grm-port-r11-20261105#$F/port#g; s#/home/layth/review/grm-boat-runner-20261031#$F/runner#g" > "$WK/rollback.sh"
newall() { # every target hashes to its next/ file
  local ok=yes
  cmp -s "$F/port/drive/drive.sh" "$HB/next/drive.sh.next" || ok=no; cmp -s "$F/runner/grm-boat.sh" "$HB/next/grm-boat.sh.next" || ok=no
  for f in mtb-box.sh r11-box.sh stage-manifest.py chain-lib.sh h1-box.sh h1-manifest.py; do cmp -s "$F/port/kit/$f" "$HB/next/kit/$f.next" || ok=no; done
  echo $ok; }
mshim() { # NAME CMD PATTERN: a shim of CMD that fails when an argument matches PATTERN (a rollback's temp file passes)
  mkdir -p "$WK/$1"; printf '#!/usr/bin/env bash\nfor a in "$@"; do case "$a" in *.kit10-back-*) ;; %s) echo "%s: refused (test)" >&2; exit 1;; esac; done\nexec %s "$@"\n' \
    "$3" "$2" "$(command -v "$2")" > "$WK/$1/$2"; chmod +x "$WK/$1/$2"; }

# I1 (review case A): a clean install
reset; inst i1; rc=$?
res I1-clean-install "0,KIT10 installed: 8 files, syntax ok,yes,0" "$rc,$(last i1),$(newall),$(temps)" ""
res I1-modes-from-the-replaced-or-a-neighbour "$(stat -c %a "$O/port/kit/r11-box.sh"),$(stat -c %a "$O/port/kit/r11-box.sh"),$(stat -c %a "$O/port/kit/stage-manifest.py"),$(stat -c %a "$O/port/drive/drive.sh")" \
  "$(stat -c %a "$F/port/kit/chain-lib.sh"),$(stat -c %a "$F/port/kit/h1-box.sh"),$(stat -c %a "$F/port/kit/h1-manifest.py"),$(stat -c %a "$F/port/drive/drive.sh")" ""
res I1-kit-push-list "chain-lib.sh ckpt.sh guard.sh h1-box.sh h1-manifest.py mtb-box.sh pg-run.sh r11-box.sh stage-manifest.py" "$(push)" ""
res I1-pre-kit10-holds-the-five-replaced "drive.sh grm-boat.sh mtb-box.sh r11-box.sh stage-manifest.py" "$(ls "$HB/pre-kit10" 2> /dev/null | LC_ALL=C sort | paste -sd' ')" ""
# I2: INSTALL.md's rollback, as written: the live files back (modes too), its hash check passes, pre-kit10/ gone
(cd "$HB" && bash "$WK/rollback.sh") > "$WK/out.i2" 2>&1; rc=$?
res I2-install-md-rollback-restores-and-checks "0,yes,no,5" "$rc,$(yn same),$(yn test -e "$HB/pre-kit10"),$(grep -c ': OK$' "$WK/out.i2")" "$(grep -v ': OK$' "$WK/out.i2" | head -n 1 | cut -c1-100)"
# I3: install again after that rollback, with nothing else done by hand
inst i3; rc=$?
res I3-install-after-the-rollback 0 "$rc" "$(last i3)"
# I4: a second run changes nothing
a=$(tree "$F"); inst i4; rc=$?
res I4-second-run-refused-tree-unchanged "1,yes" "$rc,$(yn test "$a" = "$(tree "$F")")" "$(grep -m1 STOP "$WK/out.i4" | cut -c1-100)"
# I5: a drifted live file
reset; echo "# local edit" >> "$F/runner/grm-boat.sh"; a=$(tree "$F"); inst i5; rc=$?
res I5-drifted-live-file-refused-tree-unchanged "1,yes,no" "$rc,$(yn test "$a" = "$(tree "$F")"),$(yn test -e "$HB/pre-kit10")" "$(grep -m1 STOP "$WK/out.i5" | cut -c1-100)"
# I6: a changed handback file
reset; H2=$WK/HB2; rm -rf "$H2"; cp -a "$HB" "$H2"; echo "# x" >> "$H2/next/kit/h1-box.sh.next"; inst i6 "$H2"; rc=$?
res I6-changed-handback-file-refused-tree-unchanged "1,yes" "$rc,$(yn same)" "$(grep -m1 STOP "$WK/out.i6" | cut -c1-100)"
# I7 (review case B): the runner's dir is read only. Round 1 had replaced drive.sh before the runner failed.
reset; m=$(stat -c %a "$F/runner"); chmod a-w "$F/runner"; inst i7; rc=$?; chmod "$m" "$F/runner"
# step 1 must stop it: no file installed and then rolled back (the rollback alone would also leave the tree as it was)
res I7-read-only-target-dir-stops-before-any-change "1,yes,no,0,0" "$rc,$(yn same),$(yn test -e "$HB/pre-kit10"),$(grep -c '^installed ' "$WK/out.i7"),$(grep -c '^ROLLED BACK' "$WK/out.i7")" "$(grep -m1 'STOP' "$WK/out.i7" | cut -c1-110)"
# I8: a failure in step 2 after drive.sh is in (a chmod shim refuses grm-boat.sh's temp file): all rolled back
reset; mshim sh8 chmod '*grm-boat.sh.kit10-*'; SHIM=$WK/sh8 inst i8; rc=$?
res I8-step-2-failure-rolls-back "1,yes,no,0,1" "$rc,$(yn same),$(yn test -e "$HB/pre-kit10"),$(temps),$(grep -c '^ROLLED BACK: ' "$WK/out.i8")" "$(last i8)"
inst i8b; rc=$?
res I8b-rerun-after-the-rollback 0 "$rc" "$(last i8b)"
# I9 (review case C): chmod fails on h1-box.sh's temp file: no temp file left in the kit repo's push list
reset; mshim sh9 chmod '*h1-box.sh.kit10-*'; SHIM=$WK/sh9 inst i9; rc=$?
res I9-failure-after-mktemp-leaves-no-temp "1,yes,0,ckpt.sh guard.sh mtb-box.sh pg-run.sh r11-box.sh stage-manifest.py" "$rc,$(yn same),$(temps),$(push)" "$(last i9)"
# I10: the last file's move fails: the five replaced files come back from pre-kit10/ and the two new ones go
reset; mshim sh10 mv '*h1-manifest.py.kit10-*'; SHIM=$WK/sh10 inst i10; rc=$?
res I10-last-move-fails-all-seven-undone "1,yes,no,1" "$rc,$(yn same),$(yn test -e "$HB/pre-kit10"),$(grep -c '^ROLLED BACK: ' "$WK/out.i10")" "$(last i10)"
# I11: a drive running (bash .../drive.sh) stops the install; I12: a command line that only names drive.sh doesn't
mkdir -p "$WK/run"; printf 'sleep 30\n' > "$WK/run/drive.sh"
reset; bash "$WK/run/drive.sh" & dp=$!; sleep 0.3; inst i11; rc=$?; kill "$dp" 2> /dev/null; wait "$dp" 2> /dev/null
res I11-running-drive-stops-the-install "1,yes,1" "$rc,$(yn same),$(grep -c 'a drive, runner or chain is running' "$WK/out.i11")" ""
reset; python3 -I -c 'import time; time.sleep(30)' "$WK/run/drive.sh" & pp=$!; sleep 0.3; inst i12; rc=$?; kill "$pp" 2> /dev/null; wait "$pp" 2> /dev/null
res I12-a-command-naming-drive-sh-does-not-stop-it 0 "$rc" "$(grep -m1 STOP "$WK/out.i12" | cut -c1-100)"
# I14-I17 (REVIEW-KIT10-R2 L1): a signal in step 2. A shim holds h1-box.sh's temp file for 3 s, in its chmod (I14 to
# I16) or just after its mv (I17), and drops a flag; then the signal goes in. The install has to roll back: rc 1, the
# tree as before, no pre-kit10/ and no temp file. INT goes to the install's process group, as a terminal sends it.
# Job control (set -m) gives the install its own group and keeps a background job from ignoring INT.
hshim() { # NAME CMD PATTERN pre|post: a CMD shim that, for an argument matching PATTERN, drops $WK/NAME.flag and
  # sleeps 3 s, before (pre) or after (post) it runs the real CMD
  local real; real=$(command -v "$2"); mkdir -p "$WK/$1"
  if [ "$4" = pre ]; then
    printf '#!/usr/bin/env bash\nfor a in "$@"; do case "$a" in %s) echo $$ > %q; sleep 3;; esac; done\nexec %s "$@"\n' "$3" "$WK/$1.flag" "$real"
  else
    printf '#!/usr/bin/env bash\n%s "$@"; rc=$?\nfor a in "$@"; do case "$a" in %s) echo $$ > %q; sleep 3;; esac; done\nexit $rc\n' "$real" "$3" "$WK/$1.flag"
  fi > "$WK/$1/$2"; chmod +x "$WK/$1/$2"; }
sig() { # OUT SIGNAL SHIM: the install in the background with the shim; SIGNAL once the flag is there; then wait for
  # the install and the shim, and give the install's rc
  local s=$2 d=$WK/$3 pid rc i
  rm -f "$d.flag"
  set -m
  (cd "$HB" && KIT10_D=$F/port KIT10_R=$F/runner PATH=$d:$PG:$PATH exec bash install-kit10.sh) > "$WK/out.$1" 2>&1 &
  pid=$!
  set +m
  for i in $(seq 100); do [ -s "$d.flag" ] && break; sleep 0.1; done
  if [ "$s" = INT ]; then kill -INT -- -"$pid"; else kill -"$s" "$pid"; fi
  wait "$pid" 2> /dev/null; rc=$?
  while [ -s "$d.flag" ] && kill -0 "$(cat "$d.flag")" 2> /dev/null; do sleep 0.1; done
  return $rc; }
hshim sh14 chmod '*h1-box.sh.kit10-*' pre; hshim sh17 mv '*h1-box.sh.kit10-*' post
for c in "I14 TERM sh14 term" "I15 HUP sh14 hup" "I16 INT sh14 int-to-the-group"; do
  set -- $c; reset; sig "$1" "$2" "$3"; rc=$?
  res "$1-$4-in-step-2-rolls-back" "1,yes,no,0,1" "$rc,$(yn same),$(yn test -e "$HB/pre-kit10"),$(temps),$(grep -c '^STOP: interrupted by a signal' "$WK/out.$1")" "$(last "$1")"
done
inst i16b; rc=$?
res I16b-rerun-after-an-interrupted-install 0 "$rc" "$(last i16b)"
reset; sig I17 TERM sh17; rc=$?
res I17-term-just-after-a-move-rolls-that-file-back "1,yes,no,0,1" "$rc,$(yn same),$(yn test -e "$HB/pre-kit10"),$(temps),$(grep -c '^STOP: interrupted by a signal' "$WK/out.I17")" "$(last I17)"
# I18: INSTALL.md's rollback after an install that left temp files beside two targets
reset; inst i18; cp "$F/port/kit/h1-box.sh" "$F/port/kit/h1-box.sh.kit10-AbC123"; cp "$F/runner/grm-boat.sh" "$F/runner/grm-boat.sh.kit10-XyZ789"
(cd "$HB" && bash "$WK/rollback.sh") > "$WK/out.i18" 2>&1; rc=$?
res I18-install-md-rollback-removes-left-temp-files "0,yes,0" "$rc,$(yn same),$(temps)" "$(last i18)"
reset; rm -rf "$HB/pre-kit10"
res I13-fake-tree-put-back yes "$(yn same)" ""
echo "FAILS $FAILS"
[ "$FAILS" = 0 ]
