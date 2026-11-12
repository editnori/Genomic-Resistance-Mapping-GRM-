#!/usr/bin/env python3
"""KIT11 mutants: each is one exact replacement in a copy of KIT11's runner or drive (cand/), and the KIT11 test that
covers it must name the given check among its DIFFERS. A mutant whose old text is missing or appears more than once
is an error, not a pass. Writes one line per mutant to OUT, then "MUTANTS killed/total", and exits 1 unless every
mutant is killed.
usage: python3 -I -B mutants-kit11.py HARNESS WORK OUT [JOBS]"""
import concurrent.futures
import os
import shutil
import subprocess
import sys

R, D = "cand/grm-boat.sh", "cand/drive.sh"
MUT = [
    # stage-local's box side
    ("sl-no-manifest-recheck", R, '[ "$m" = "$want" ] || bye "the manifest copy is not the pinned one" 3', "true",
     "T16-manifest-changed-on-the-way-86"),
    ("sl-source-link-allowed", R, '[ -d "$s" ] && [ ! -L "$s" ] || bye', '[ -d "$s" ] || bye', "T11-source-folder-a-link-86"),
    ("sl-follows-folder-links", R, "F = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW", "F = os.O_RDONLY | os.O_DIRECTORY",
     "T10-folder-link-in-source-86"),
    ("sl-no-regular-file-check", R, '            if not stat.S_ISREG(st.st_mode):\n                print("LINKERR not a regular',
     '            if False:\n                print("LINKERR not a regular', "T9-names-the-link"),
    ("sl-copies-not-links", R, "os.link(parts[-1], os.path.join(dst, p), src_dir_fd=fd, follow_symlinks=False)",
     'open(os.path.join(dst, p), "wb").write(open("/proc/self/fd/%d/%s" % (fd, parts[-1]), "rb").read())',
     "T1-hard-links-not-copies"),
    ("sl-dest-link-allowed", R, '[ ! -L "$d" ] || bye "the destination is a link" 4', "true", "T15-dest-a-link-86"),
    ("sl-keeps-dest-on-failure", R, 'made=1\npython3 -I - "$s" "$d" "$t"', 'made=0\npython3 -I - "$s" "$d" "$t"',
     "T12-changed-bytes-86-source-kept"),
    ("sl-failure-removes-the-source", R,
     'bye() { rm -f "$t"; [ $made = 1 ] && rm -rf "$d"; echo "REFUSED $1"; exit $2; }\nm=$(sha256sum "$t" | cut -d" " -f1); echo "TSHA $m"',
     'bye() { rm -f "$t"; [ $made = 1 ] && rm -rf "$d" "$s"; echo "REFUSED $1"; exit $2; }\nm=$(sha256sum "$t" | cut -d" " -f1); echo "TSHA $m"',
     "T12-changed-bytes-86-source-kept"),
    # stage-local's checks here
    ("pd-check-only-takes-receipts", R, "        if is_receipt(p) or combine.fullmatch(p):", "        if False:",
     "T4-receipt-goes-through-push-data-97"),
    ("pd-check-only-fills-combine", R, '    if combine_dest: problems.append("stage-local does not fill /home/user/grm-combine")\n', "",
     "T5-grm-combine-refused-97"),
    ("sl-dest-inside-source", R, '  case "$dest/" in "$src"/*) fail 2 "$ID" stage-local "--dest $dest is --from $src or inside it";; esac\n', "",
     "T7-folders-must-be-plain-and-apart"),
    ("sl-source-inside-dest", R, '  case "$src/" in "$dest"/*) fail 2 "$ID" stage-local "--from $src is inside --dest $dest";; esac\n', "",
     "T7-folders-must-be-plain-and-apart"),
    ("sl-tmp-allowed", R, '[[ "$x" =~ ^/home/user(/[A-Za-z0-9._-]+)+$ ]]', '[[ "$x" =~ ^/(home/user|tmp)(/[A-Za-z0-9._-]+)+$ ]]',
     "T7-folders-must-be-plain-and-apart"),
    ("sl-no-pin-check-here", R, '  out=$(pd_tool "$src" "$STAGE/manifest-$ts.sha256" "$want" - "$dest" 2>&1); rc=$?',
     '  out=$(pd_tool "$src" "$STAGE/manifest-$ts.sha256" "$(sha256sum "$STAGE/manifest-$ts.sha256" | cut -d" " -f1)" - "$dest" 2>&1); rc=$?',
     "T2-wrong-pin-97-no-box-call"),
    ("pd-box-skips-the-shared-checks", R, "'\"$PD_CHECK\"\ncmd_push_data() {", "'\ncmd_push_data() {",
     "L6-linux-push-data-stages-here"),
    # the Linux host mode
    ("lin-no-floor", R, '    [ "$free" -ge "$MIN_AVAIL_MB" ] && [ "$commit" -ge "$MIN_AVAIL_MB" ] ||', "    true ||",
     "L4-linux-memory-under-the-floor-76"),
    ("lin-stages-through-wslpath", R, 'if [ "${HOST_KIND:-wsl}" = linux ]; then base=$STAGE_DIR/stage-$LANE; else',
     "if false; then base=$STAGE_DIR/stage-$LANE; else", "L6-linux-push-data-stages-here"),
    ("lin-no-boat-needed", R,
     '  linux) [ -n "${GRM_BOAT:-}" ] || { echo "grm-boat: GRM_HOST=linux needs GRM_BOAT, the boat binary" >&2; exit 2; };;',
     "  linux) ;;", "L1-linux-needs-grm-boat"),
    ("lin-any-host", R, '  *) echo "grm-boat: GRM_HOST must be wsl or linux, not $HOST_KIND" >&2; exit 2;;', "  *) ;;",
     "L2-unknown-host-2"),
    ("lin-reads-windows-memory", R, '  if [ "$HOST_KIND" = linux ]; then\n    [[ "$MIN_AVAIL_MB"', '  if false; then\n    [[ "$MIN_AVAIL_MB"',
     "L3-linux-preflight-reads-this-host"),
    # drive.sh's snapshot mode
    ("drv-always-push-data", D, '    if [ $SNAPDATA = 1 ]; then waited "$R" stage-local', '    if false; then waited "$R" stage-local',
     "D1-three-stage-local-no-push-data"),
    ("drv-snapshot-any-kind", D, '  [ "$KIND" = h1 ] || { echo "drive: DATA_SNAPSHOT_DIR is for KIND=h1 only" >&2; exit 2; }', "  :",
     "D3-snapshot-mode-is-h1-only"),
    ("drv-snapshot-no-template", D,
     '  case "${TEMPLATE:-}" in grm-tmpl-*) ;; *) echo "drive: DATA_SNAPSHOT_DIR needs TEMPLATE, the grm-tmpl-* data snapshot" >&2; exit 2;; esac',
     "  :", "D4-snapshot-mode-needs-a-grm-tmpl-template"),
    ("drv-snapshot-vocab-from-root", D, "    if [ $SNAPDATA = 1 ]; then\n      # KIT11: the files are on the snapshot",
     "    if false; then\n      # KIT11: the files are on the snapshot", "D1-snapshot-pilot-rc"),
    ("drv-snapshot-vocab-unpinned", D, '[ "$VSHA" = "$VOCAB_MANIFEST_SHA256" ] || die 2 "the vocabulary manifest from',
     'true || die 2 "the vocabulary manifest from', "D6-other-vocab-manifest-2-no-box"),
    ("drv-snapshot-folder-anywhere", D, '[[ "${!v:-}" =~ ^/home/user(/[A-Za-z0-9._-]+)+$ ]] || { echo "drive: $v must be',
     'true || { echo "drive: $v must be', "D7-snapshot-folders-under-home-user"),
    ("drv-snapshot-needs-data-root", D, '[ $SNAPDATA = 1 ] || need="$need DATA_ROOT"', 'need="$need DATA_ROOT"',
     "D10-snapshot-mode-needs-no-data-or-vocab-root"),
    ("drv-push-mode-no-data-root", D, '[ $SNAPDATA = 1 ] || need="$need DATA_ROOT"', ":", "D8-push-mode-still-needs-data-root"),
    ("drv-snapshot-no-vocab-manifest", D,
     '  [ -f "${VOCAB_MANIFEST:-}" ] || { echo "drive: DATA_SNAPSHOT_DIR needs VOCAB_MANIFEST, the vocabulary manifest file" >&2; exit 2; }',
     "  :", "D5-missing-vocab-manifest-2-no-box"),
]


def run(h, work, m):
    mid, f, old, new, want = m
    d = os.path.join(work, mid)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    shutil.copytree(os.path.join(h, "cand"), os.path.join(d, "cand"))
    p = os.path.join(d, f)
    text = open(p).read()
    if text.count(old) != 1:
        return "%s ERROR old text found %d times in %s" % (mid, text.count(old), f)
    open(p, "w").write(text.replace(old, new))
    if f == R:
        cmd = ["bash", os.path.join(h, "tests/test-runner-kit11.sh"), p, os.path.join(d, "w")]
    else:
        cmd = ["bash", os.path.join(h, "tests/test-drive-h1-kit11.sh"), os.path.join(d, "w"), p, os.path.join(h, "kit")]
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
