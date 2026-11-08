#!/usr/bin/env python3
"""KIT12 mutants: each is one exact replacement in a copy of KIT12's drive (cand/drive.sh) or kit (kit/h1-box.sh,
kit/r11-box.sh), and the KIT12 test that covers it must name the given check among its DIFFERS. A mutant whose old
text is missing or appears more than once is an error, not a pass. Writes one line per mutant to OUT, then
"MUTANTS killed/total", and exits 1 unless every mutant is killed.
usage: python3 -I -B mutants-kit12.py HARNESS WORK OUT [JOBS]"""
import concurrent.futures
import os
import shutil
import subprocess
import sys

D, HB, RB = "cand/drive.sh", "kit/h1-box.sh", "kit/r11-box.sh"
MUT = [
    # drive.sh's addon role
    ("drv-addon-is-a-run", D, "addon) H1ROLE=addon;;", "addon) H1ROLE=run;;", "A2-a-run-box-s-data-then-the-main-part"),
    ("drv-any-part-cell", D, '[[ "${PARTS:-}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || die 2 "stage $BOX: an addon row',
     'true || die 2 "stage $BOX: an addon row', "A6-each-names-its-cell"),
    ("drv-main-any-kind", D, "      pilot|s1|s2|conf) ;;\n", "      pilot|s1|s2|conf|addon|combine) ;;\n",
     "A6-bad-main-part-cells-2-no-box"),
    ("drv-addon-no-binary-order", D, 'if [ "$KIND" = h1 ] && [ "$H1ROLE" != combine ]; then\n  [ -n "$REF_BOX" ]',
     'if [ "$KIND" = h1 ] && [ "$H1ROLE" = run ]; then\n  [ -n "$REF_BOX" ]', "A9-addon-waits-for-the-reference-binary"),
    ("drv-main-not-done-ok", D, '[ -e "$MP/DONE" ] && [ -d "$MP/out/h1work" ] || die 79', '[ -d "$MP/out/h1work" ] || die 79',
     "A13-main-part-pulled-but-not-done-79"),
    ("drv-main-not-pulled-ok", D, '[ -e "$MP/DONE" ] && [ -d "$MP/out/h1work" ] || die 79', '[ -e "$MP/DONE" ] || die 79',
     "A12-main-part-done-but-not-pulled-79"),
    ("drv-flags-unchecked", D, '      cmp -s "$mf" "$SD/H1-FLAGS" ||\n', "      true ||\n", "A7-other-flags-than-the-main-part-2-no-box"),
    ("drv-flags-unread", D, '      { [ -f "$mf" ] && [ "$(wc -l < "$mf")" = 1 ] && grep -qxE', '      { true || grep -qxE',
     "A8-main-part-without-flags-2-no-box"),
    ("drv-set-not-fresh", D, 'rm -rf "$SD/ref-in"; H1REF=', "H1REF=", "A14-rerun-builds-the-set-afresh"),
    ("drv-set-unpinned-in-notes", D, "      say \"main part $PARTS's records pinned:", "      : \"main part $PARTS's records pinned:",
     "A2-pinned-in-notes-before-the-box"),
    ("drv-no-records-push", D, '    if [ "$H1ROLE" = addon ]; then\n      # KIT12: run records go through push-data',
     '    if false; then\n      # KIT12: run records go through push-data', "A2-a-run-box-s-data-then-the-main-part"),
    ("drv-records-other-pin", D, '--manifest "$SD/REF-MANIFEST.sha256" --manifest-sha256 "$RSHA"',
     '--manifest "$SD/REF-MANIFEST.sha256" --manifest-sha256 "$MSHA"', "A2-the-main-part-s-records-are-sent"),
    ("drv-records-refusal-ignored", D, "        bail $? \"push-data of main part $PARTS's records refused or failed",
     "        : \"push-data of main part $PARTS's records refused or failed", "A10-records-refused-97-closed"),
    ("drv-env-no-reference", D, "H1_REFERENCE=${H1REF:-}\n", "H1_REFERENCE=\n", "A2-stage-env-names-the-part-and-pin"),
    ("drv-env-no-reference-pin", D, "H1_REFERENCE_SHA256=${RSHA:-}", "H1_REFERENCE_SHA256=", "A2-stage-env-names-the-part-and-pin"),
    ("drv-env-no-meta-pin", D, 'H1_META_SHA256=$([ "$H1ROLE" != combine ] && echo "$MSHA")',
     'H1_META_SHA256=$([ "$H1ROLE" = run ] && echo "$MSHA")', "A2-stage-env-names-the-part-and-pin"),
    ("drv-addon-inputs-as-combine", D, '  if [ "$H1ROLE" != combine ]; then\n    python3 -I "$KIT/h1-manifest.py" meta',
     '  if [ "$H1ROLE" = run ]; then\n    python3 -I "$KIT/h1-manifest.py" meta', "A2-pinned-in-notes-before-the-box"),
    # h1-box.sh's addon role
    ("hb-union-no-reference", HB, '"${EXTRA[@]}" "${REFARG[@]}" --estimate-only --free-disk-bytes 1000000000000000',
     '"${EXTRA[@]}" --estimate-only --free-disk-bytes 1000000000000000', "B1-addon-union-passes-the-reference"),
    ("hb-verify-no-reference", HB, '"${EXTRA[@]}" "${REFARG[@]}" --estimate-only --free-disk-bytes "$use"',
     '"${EXTRA[@]}" --estimate-only --free-disk-bytes "$use"', "B2-addon-verify-passes-the-reference"),
    ("hb-run-no-reference", HB, '    "${TRACE[@]}" "${EXTRA[@]}" "${REFARG[@]}"\n', '    "${TRACE[@]}" "${EXTRA[@]}"\n',
     "B3-addon-run-passes-the-reference"),
    ("hb-union-unchecked-set", HB, 'checkset "$VOC" "$H1_VOCAB_SHA256"; refset\n', 'checkset "$VOC" "$H1_VOCAB_SHA256"\n',
     "B5-reference-pin-mismatch-2-no-run"),
    ("hb-run-unchecked-set", HB, 'RC=2 fail "RAW is missing or writable; run verify first"\n  refset\n',
     'RC=2 fail "RAW is missing or writable; run verify first"\n', "B12-run-rechecks-the-set"),
    ("hb-other-parts-ok", HB, '-printf \'%f\\n\' 2> /dev/null)" = "$H1_REFERENCE" ] ||', '-printf \'%f\\n\' 2> /dev/null)" ] ||',
     "B7-a-second-part-2"),
    ("hb-unfinished-ok", HB, '[ -d "$REFW" ] && [ -f "$REFW/summary.json" ] ||', '[ -d "$REFW" ] ||',
     "B10-unfinished-main-part-2"),
    ("hb-any-reference-name", HB, '[[ "${H1_REFERENCE:-}" =~ ^[a-z0-9][a-z0-9._-]*$ ]] ||', '[[ -n "${H1_REFERENCE:-}" ]] ||',
     "B9-reference-name-a-path-2"),
    ("hb-no-addon-role", HB, 'if [ "${H1_ROLE:-run}" = addon ]; then\n  [[ "${H1_REFERENCE', 'if false; then\n  [[ "${H1_REFERENCE',
     "B1-addon-union-passes-the-reference"),
    ("hb-reference-on-any-role", HB, 'if [ "${H1_ROLE:-run}" = addon ]; then\n  [[ "${H1_REFERENCE',
     'if [ -n "${H1_REFERENCE:-}" ]; then\n  [[ "${H1_REFERENCE', "B11-a-run-box-passes-no-reference"),
    # r11-box.sh's dispatch
    ("rb-addon-command", RB, 'H1CMD=${H1_ROLE:-run}; [ "$H1CMD" = addon ] && H1CMD=run', "H1CMD=${H1_ROLE:-run}",
     "D2-addon-start-runs-h1-box-run-to-done"),
    ("rb-addon-pause-file", RB, 'gp=$(tail -n 1 "$W/gates/h1-$H1CMD.tripped"', 'gp=$(tail -n 1 "$W/gates/h1-${H1_ROLE:-run}.tripped"',
     "D3-addon-reads-h1-run-s-pause"),
    ("rb-addon-not-a-role", RB, "run|addon|combine) ;;", "run|combine) ;;", "D2-addon-start-runs-h1-box-run-to-done"),
]


def run(h, work, m):
    mid, f, old, new, want = m
    d = os.path.join(work, mid)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    shutil.copytree(os.path.join(h, "cand"), os.path.join(d, "cand"))
    shutil.copytree(os.path.join(h, "kit"), os.path.join(d, "kit"), symlinks=True)
    p = os.path.join(d, f)
    text = open(p).read()
    if text.count(old) != 1:
        return "%s ERROR old text found %d times in %s" % (mid, text.count(old), f)
    open(p, "w").write(text.replace(old, new))
    if f == D:
        cmd = ["bash", os.path.join(h, "tests/test-drive-h1-kit12.sh"), os.path.join(d, "w"), p, os.path.join(d, "kit")]
    else:
        cmd = ["bash", os.path.join(h, "tests/test-h1-box-kit12.sh"), os.path.join(d, "w"), os.path.join(d, "kit")]
    o = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=1800).stdout.decode("utf-8", "replace")
    open(os.path.join(d, "test.out"), "w").write(o)
    differs = [x.split()[1] for x in o.splitlines() if x.startswith("RESULT ") and " DIFFERS" in x]
    verdict = "KILLED" if want in differs else "SURVIVED"
    return "%s %s by %s (DIFFERS: %s)" % (mid, verdict, want, " ".join(differs) or "none")


def main():
    h, work, out = sys.argv[1:4]
    jobs = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    ids = [m[0] for m in MUT]
    assert len(ids) == len(set(ids)), "a mutant id twice"
    with concurrent.futures.ThreadPoolExecutor(jobs) as ex:
        lines = list(ex.map(lambda m: run(h, work, m), MUT))
    killed = sum(1 for x in lines if " KILLED " in x)
    with open(out, "w") as f:
        for x in lines:
            f.write(x + "\n")
        f.write("MUTANTS %d/%d\n" % (killed, len(MUT)))
    print(open(out).read(), end="")
    sys.exit(0 if killed == len(MUT) else 1)


main()
