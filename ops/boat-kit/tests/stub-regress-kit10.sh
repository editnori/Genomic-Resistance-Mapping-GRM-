#!/usr/bin/env bash
# KIT10 (from KIT9's stub-regress.sh): PORT's 32 stub cases (drv-gates-v4.sh: stub runner, stub boat, no box) on the
# live drive.sh v7 and the candidate; only RESULT lines are compared. One harness change, for both drives: PORT's stub
# runner gains an extend that prints "END <now + ttl>" and exits 0, as grm-boat.sh extend does (KIT10 round 2,
# REVIEW-KIT10 L2: the candidate takes its deadline from that line and stops with 86 without it). The candidate sets a
# slow host's lifetime through the runner's extend where v7 called boat extend directly; the stub boat already
# answered extend with 0. The candidate is also run on the unchanged harness, and its differences are printed
# (expected: the slow-host case only).
# Exit 1 unless both give 32 cases, the same results, and 30 as-predicted plus the 2 known K2 DIFFERS.
# usage: stub-regress-kit10.sh WORK CANDIDATE_DRIVE_SH
set -u
W=${1:?work} NEW=${2:?candidate drive.sh}
# KIT10: drv-gates' D1b and D4 cases send INT, TERM and HUP to drive.sh, which can't trap a signal ignored when it
# started, and D4 waits for it without a timeout. A chain started with nohup or a background & from a script ignores
# HUP or INT (run 3 waited on D4-INT from 12:07Z), so refuse at once instead.
ign=$(sed -n 's/^SigIgn:[[:space:]]*//p' /proc/$$/status 2> /dev/null)
[ -n "$ign" ] && (( 0x$ign & 0x4003 )) && { echo "STUB-REGRESS FAIL: HUP, INT or TERM is ignored here (SigIgn $ign); start the chain with setsid -f, without nohup or &"; exit 1; }
PD=/home/layth/review/grm-port-r11-20261105 V7=/home/layth/review/grm-kit9-20261106/live/drive.sh
rm -rf "$W"; mkdir -p "$W"
# the extend line goes in before the stub's catch-all; after the stub's shift its arguments are ID --ttl N
EXT='  extend) echo "END $(( $(date +%s) + $3 ))"; exit 0;;'
python3 -I - "$PD/dryrun/drv-stub-runner-v4.sh" "$W/drv-stub-runner-ext.sh" "$EXT" << 'EOF'
import sys
src, out, ext = sys.argv[1:]
t = open(src).read(); a = '  *) echo "stub runner: $c?" >&2; exit 2;;\n'
if t.count(a) != 1: sys.exit("the stub runner has %d catch-all lines, not 1" % t.count(a))
open(out, "w").write(t.replace(a, ext + "\n" + a))
EOF
chmod +x "$W/drv-stub-runner-ext.sh"; grep -qxF "$EXT" "$W/drv-stub-runner-ext.sh" || { echo "STUB-REGRESS FAIL: no extend line"; exit 1; }
sed "s#\$ME/drv-stub-runner-v4.sh#$W/drv-stub-runner-ext.sh#" "$PD/dryrun/drv-gates-v4.sh" > "$W/drv-gates-ext.sh"
grep -q "RUNNER=$W/drv-stub-runner-ext.sh" "$W/drv-gates-ext.sh" || { echo "STUB-REGRESS FAIL: harness not pointed at the extended stub"; exit 1; }
sha256sum "$V7" "$NEW" "$PD/dryrun/drv-gates-v4.sh" "$PD/dryrun/drv-stub-runner-v4.sh" "$W/drv-gates-ext.sh" "$W/drv-stub-runner-ext.sh"
diff "$PD/dryrun/drv-stub-runner-v4.sh" "$W/drv-stub-runner-ext.sh"
bash "$W/drv-gates-ext.sh" "$W/v7" "$V7" > "$W/v7-drv-gates.out" 2>&1
bash "$W/drv-gates-ext.sh" "$W/new" "$NEW" > "$W/new-drv-gates.out" 2>&1
bash "$PD/dryrun/drv-gates-v4.sh" "$W/new-plain" "$NEW" > "$W/new-plain-drv-gates.out" 2>&1
r() { grep -o '^RESULT [^ ]* want=[^ ]* got=[^ ]* [^ ]*' "$1"; }
a=$(r "$W/v7-drv-gates.out" | wc -l) b=$(r "$W/new-drv-gates.out" | wc -l)
echo "drv-gates (extended stub): v7 $a cases, candidate $b cases; differences:"
diff <(r "$W/v7-drv-gates.out") <(r "$W/new-drv-gates.out"); d=$?
[ $d = 0 ] && echo none
echo "candidate on the unchanged harness, against v7 (information):"
diff <(r "$W/v7-drv-gates.out") <(r "$W/new-plain-drv-gates.out") && echo none
echo "candidate results:"; r "$W/new-drv-gates.out" | awk '{print $NF}' | sort | uniq -c
ap=$(r "$W/new-drv-gates.out" | grep -c ' as-predicted$')
[ "$a" = 32 ] && [ "$b" = 32 ] && [ $d = 0 ] && [ "$ap" = 30 ] && { echo "STUB-REGRESS PASS"; exit 0; }
echo "STUB-REGRESS FAIL"; exit 1
