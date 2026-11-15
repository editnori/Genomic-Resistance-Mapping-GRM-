#!/usr/bin/env python3
"""KIT13's planted faults: six wrong versions of drive.sh's two-slow-hosts stop, each run under
test-drive-h1-kit13.sh on the stub runner and stub boat (no box). Each fault is one exact text replacement, which must
match once. A fault is caught when the test fails, and the DIFFERS ids say which checks caught it. The unchanged
drive.sh must pass first.

usage: python3 -I -B mutants-kit13.py DRIVE_SH KIT_DIR TEST_SH WORK OUT_JSON
"""
import json
import pathlib
import re
import subprocess
import sys

STOP = '    die 78 "budget: the last two boxes had slow hosts; stop and report"'
PASS_IF = '  if [ "${SLOW_PAIR_PASS:-0}" = 1 ] && [ "$SLOW_OK" = 1 ]; then'
NOTE = ('    say "budget: the last two boxes had slow hosts; SLOW_PAIR_PASS=1 and stage $BOX is SLOW_OK, '
        'so this drive goes on"')
CHECK = '  *) echo "drive: SLOW_PAIR_PASS is 0 or 1, not ${SLOW_PAIR_PASS}" >&2; exit 2;;'
CAP = ("awk -v u=\"$used\" -v e=\"$EST_SLOW\" -v c=\"$BUDGET_CAP\" 'BEGIN { exit !(u + e > c) }' && die 78")

FAULTS = [
    ("F1-pass-ignores-slow-ok", PASS_IF, '  if [ "${SLOW_PAIR_PASS:-0}" = 1 ]; then'),
    ("F2-pass-is-the-default", PASS_IF, '  if [ "${SLOW_PAIR_PASS:-1}" = 1 ] && [ "$SLOW_OK" = 1 ]; then'),
    ("F3-any-value-accepted", CHECK, '  *) ;;'),
    ("F4-no-note", NOTE, '    :'),
    ("F5-pass-skips-the-cap", CAP, CAP.replace("&& die 78", '&& [ "${SLOW_PAIR_PASS:-0}" != 1 ] && die 78')),
    ("F6-stop-text-changed", STOP, STOP.replace("stop and report", "stopped")),
]


def run(drive, kit, test, work):
    p = subprocess.run(["bash", test, str(work), str(drive), kit], capture_output=True, text=True, timeout=900)
    differs = sorted(re.findall(r"^RESULT (\S+) .* DIFFERS", p.stdout, re.M))
    return p.returncode, differs


def main(argv):
    drive, kit, test, work, out = argv[1:6]
    work = pathlib.Path(work)
    work.mkdir(parents=True, exist_ok=True)
    text = pathlib.Path(drive).read_text()
    rc, differs = run(drive, kit, test, work / "unchanged")
    result = {"unchanged": {"rc": rc, "differs": differs}, "faults": []}
    if rc != 0:
        print(f"the unchanged drive fails: {differs}")
    caught = 0
    for name, old, new in FAULTS:
        if text.count(old) != 1:
            sys.exit(f"{name}: the text to replace occurs {text.count(old)} times, not once")
        path = work / f"{name}.sh"
        path.write_text(text.replace(old, new))
        frc, fdiffers = run(path, kit, test, work / name)
        result["faults"].append({"fault": name, "rc": frc, "caught_by": fdiffers})
        caught += frc != 0 and bool(fdiffers)
        print(f"{name}: rc={frc} caught by {', '.join(fdiffers) or 'nothing'}")
    result["caught"] = f"{caught} of {len(FAULTS)}"
    pathlib.Path(out).write_text(json.dumps(result, indent=1) + "\n")
    print(f"caught {caught} of {len(FAULTS)}; unchanged rc={rc}")
    return 0 if rc == 0 and caught == len(FAULTS) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
