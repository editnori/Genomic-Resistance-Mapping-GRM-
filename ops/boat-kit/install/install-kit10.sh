#!/usr/bin/env bash
# KIT10 install on the laptop: KIT9's and KIT10's changes together (KIT10's drive.sh and grm-boat.sh hold KIT9's).
# Run it from the dir that holds next/ and SHA256SUMS. It checks every file first and changes nothing unless all
# checks pass; then each file goes through a temp file beside its target, takes the mode of the file it replaces
# (a new kit file takes a neighbour's), and is moved into place, so a script already running keeps its old file.
# Each file it replaces is kept in pre-kit10/ here, outside the kit repo, since the runner pushes untracked kit files
# to every box. KIT10_D and KIT10_R override the two roots (the tests use them).
# KIT10 round 2 (REVIEW-KIT10 L4): step 1 also checks that each target's dir is writable, so a read-only dir stops
# the install before any change. Any failure in step 2 rolls back: the temp file goes, every file already replaced
# comes back from pre-kit10/ through a temp file and mv, every new file already in goes, and the live hashes are
# checked; pre-kit10/ is removed only when all of them match. The running-script check matches only bash running
# drive.sh, grm-boat.sh or chain-tests.sh, not any command line that names one.
# usage: bash install-kit10.sh
set -u
D=${KIT10_D:-/home/layth/review/grm-port-r11-20261105} R=${KIT10_R:-/home/layth/review/grm-boat-runner-20261031}
K=$D/kit P=$PWD/pre-kit10
# new file | target | sha256 of the copy the tests ran against ("new": no file there yet) | new sha256 | mode from
L="next/drive.sh.next|$D/drive/drive.sh|db62192168156fc773c46fe101c188e3cfa2192312cb46d356a72a81f2d860f7|14ef9a19fa7af2f0a5e0ce00998a2b8b1d377f1ccce4511ab20deb47a77e27d2|$D/drive/drive.sh
next/grm-boat.sh.next|$R/grm-boat.sh|45fb0c9ff27c7a851f13c5bcbfaa1c99388a47d61f7511ce409ba69ae05c3cd2|dd3f4842c07fdf9ccc3a5da4f840066e3a1517aaea1403e4fe3ed752e9b2db5e|$R/grm-boat.sh
next/kit/mtb-box.sh.next|$K/mtb-box.sh|0a42c769c633435f7e1200522bfd2aaec53b61843935d225d56f0a8b4f850361|f67db557c7e642efb01e4a01184cc1103424e97a16e787d2d4a997569ee5f072|$K/mtb-box.sh
next/kit/r11-box.sh.next|$K/r11-box.sh|4920a98e6066fe5b39778c293ebb43f74a7a74244065e1be6795e474100f8d56|b618def77fce3d39f3d80f8ea9a6fe359850f5c77d4cdec330ff5e304b85bef1|$K/r11-box.sh
next/kit/stage-manifest.py.next|$K/stage-manifest.py|36a9be604d0cc23c515f498a46ff5d1492ea5c6b8f01774cd842d542daf99699|68cf7b4cd25fa4ddea2df407ddf51b118dda2a3dba6d2c8f5180d62c5e1eddc5|$K/stage-manifest.py
next/kit/chain-lib.sh.next|$K/chain-lib.sh|new|1e9514d0d3fdf705d043b3edb037bc889355871a76d2c96694601f1aafc28e3b|$K/r11-box.sh
next/kit/h1-box.sh.next|$K/h1-box.sh|new|b6f32caf66bef5e12403e242ed5e2666da78828f4639212919be4b8afa98e268|$K/r11-box.sh
next/kit/h1-manifest.py.next|$K/h1-manifest.py|new|12f93492cb482432960fd818e4cd843cf5d5f9ef86c4bba93326a1c6abe46e09|$K/stage-manifest.py"
h() { sha256sum < "$1" | cut -d' ' -f1; }
syntax() { # FILE...: bash -n on the shell files, an ast parse on the Python ones
  local f n=0 py=()
  for f in "$@"; do case "$f" in *.py|*.py.next) py+=("$f");; *) bash -n "$f" || { echo "SYNTAX: $f"; n=1; };; esac; done
  python3 -I -B -c 'import ast, sys
for f in sys.argv[1:]: ast.parse(open(f).read())' "${py[@]}" < /dev/null || n=1
  return $n
}

# 1. every check before any change
sha256sum -c --strict --quiet SHA256SUMS || { echo "STOP: this handback dir doesn't match its SHA256SUMS"; exit 1; }
if pgrep -af '^(/[^ ]*/)?bash( -[A-Za-z]+)* ([^ ]*/)?(drive|grm-boat|chain-tests)\.sh( |$)'; then
  echo "STOP: a drive, runner or chain is running (above); install when none is"; exit 1
fi
bad=0
while IFS='|' read -r src dst old new ref; do
  [ "$(h "$src")" = "$new" ] || { echo "STOP: $src is not $new"; bad=1; }
  if [ "$old" = new ]; then
    [ ! -e "$dst" ] || { echo "STOP: $dst already exists"; bad=1; }
  else
    [ -f "$dst" ] && [ "$(h "$dst")" = "$old" ] || { echo "STOP: $dst is not the copy the tests ran against ($old)"; bad=1; }
  fi
  [ -f "$ref" ] || { echo "STOP: no $ref to take the mode from"; bad=1; }
  [ -d "${dst%/*}" ] && [ -w "${dst%/*}" ] || { echo "STOP: can't write ${dst%/*}, so $dst can't be replaced"; bad=1; }
done <<< "$L"
[ ! -e "$P" ] || { echo "STOP: $P is already there (an earlier install?)"; bad=1; }
[ -w "$PWD" ] || { echo "STOP: can't write $PWD, so pre-kit10/ can't be made"; bad=1; }
for f in chain-lib.sh h1-box.sh h1-manifest.py; do
  git -C "$K" check-ignore -q "$f" < /dev/null && { echo "STOP: the kit repo ignores $f, so a push to a box would leave it out"; bad=1; }
done
syntax $(cut -d'|' -f1 <<< "$L") || { echo "STOP: a new file fails its syntax check (above)"; bad=1; }
[ $bad = 0 ] || exit 1

# 2. install; any failure from here rolls back what step 2 changed
t="" IN=()   # the temp file in flight; the entries whose target step 2 has replaced or made
rollback() {
  local i src dst old new ref b r ok=1
  [ -n "$t" ] && rm -f "$t"
  for (( i = ${#IN[@]} - 1; i >= 0; i-- )); do
    IFS='|' read -r src dst old new ref <<< "${IN[$i]}"
    if [ "$old" = new ]; then rm -f "$dst" || ok=0; continue; fi
    b=$P/${dst##*/}
    r=$(mktemp "$dst.kit10-back-XXXXXX") && cat "$b" > "$r" && chmod --reference="$b" "$r" && mv "$r" "$dst" ||
      { [ -n "${r:-}" ] && rm -f "$r"; echo "ROLLBACK: could not put $b back as $dst"; ok=0; }
  done
  while IFS='|' read -r src dst old new ref; do
    if [ "$old" = new ]; then [ ! -e "$dst" ] || { echo "ROLLBACK: $dst is still there"; ok=0; }
    else [ "$(h "$dst")" = "$old" ] || { echo "ROLLBACK: $dst is not $old"; ok=0; }; fi
    ls "$dst".kit10-* > /dev/null 2>&1 && { echo "ROLLBACK: a $dst.kit10-* temp file is left"; ok=0; }
  done <<< "$L"
  if [ $ok = 1 ]; then
    rm -rf "$P"
    echo "ROLLED BACK: every live file hashes as before the install and no temp file is left; pre-kit10/ removed, so a rerun starts clean"
  else
    echo "ROLLBACK INCOMPLETE (above): pre-kit10/ is kept; put the live files back by hand with INSTALL.md's rollback"
  fi
  exit 1
}
stop2() { echo "STOP: $*"; rollback; }
mkdir "$P" || { echo "STOP: could not make $P"; exit 1; }
while IFS='|' read -r src dst old new ref; do
  if [ "$old" != new ]; then cp -p "$dst" "$P/${dst##*/}" || stop2 "could not keep $dst"; fi
  t=$(mktemp "$dst.kit10-XXXXXX") || { t=""; stop2 "could not make a temp file beside $dst"; }
  cat "$src" > "$t" && chmod --reference="$ref" "$t" || stop2 "could not write $t"
  mv "$t" "$dst" || stop2 "could not move $t to $dst"
  t=""; IN+=("$src|$dst|$old|$new|$ref")
  [ "$(h "$dst")" = "$new" ] || stop2 "$dst doesn't hash to $new after the move"
  echo "installed $dst $new $(stat -c %A "$dst")"
done <<< "$L"
syntax $(cut -d'|' -f2 <<< "$L") || stop2 "an installed file fails its syntax check (above)"
echo "KIT10 installed: 8 files, syntax ok"
