#!/usr/bin/env bash
# Boot hook (run from /root/onstart.sh): bring back the pair state/active
# records on the live port, through the deployment toolkit. restore.sh with no
# pair restores state/active: it verifies the release and the publication
# against their manifests, the approval and state/active itself, resolves an
# interrupted start, starts the pair as a supervised start, requires it to own
# the port and pass smoke.sh, and records it. When that pair already serves and
# is recorded, it restarts nothing.
#
# Retry (main's ruling, 2026-10-30): a failed restore of the active pair is
# tried once more, after BOOT_RETRY_WAIT seconds (30) on top of restore.sh's
# own readiness wait, so a transient failure does not move the site to the
# fallback, which restore.sh then records in state/active for later boots.
#
# Fallback (main's ruling, 2026-10-29; moved to revision 8, 2026-10-30): when
# both tries fail, restore.sh ui-r8-690815a rev8-273-20261029 --apply, through
# the same checks, never
# serve-live.sh. Both outcomes are logged (FAIL active=..., FALLBACK ... rc=...)
# and recorded in BOOT-FALLBACK in the root. If the fallback fails too, nothing
# serves.
# usage: boot-active.sh [--dry-run]   (dry run: restore.sh without --apply,
#   the retry without its wait, and BOOT-FALLBACK is not written)
# BOOT_DEPLOY_ROOT and BOOT_SCRIPTS replace the live root and scripts (tests).
# Exit 0 when the active pair serves, 3 when the fallback serves, else 1.
set -u
export GRM_DEPLOY_ROOT=${BOOT_DEPLOY_ROOT:-/root/grm-deploy-20261025}
G=${BOOT_SCRIPTS:-/root/grm-ops-20261029/scripts/deployment}
FB_LABEL=ui-r8-690815a FB_PUB=rev8-273-20261029
BOOT_RETRY_WAIT=${BOOT_RETRY_WAIT:-30}
dry=0
[ "${1:-}" != --dry-run ] || dry=1
apply=(--apply)
((dry == 0)) || apply=()
say() { printf '%s boot-active: %s\n' "$(date -u +%FT%TZ)" "$*"; }
mark() {
    local m=$GRM_DEPLOY_ROOT/BOOT-FALLBACK
    if ((dry)); then say "dry run: would append to $m: $*"; return 0; fi
    [ ! -L "$m" ] && { [ ! -e "$m" ] || [ -f "$m" ]; } || { say "FAIL: $m is not a plain file; not written"; return 0; }
    printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >>"$m" || say "FAIL: cannot write $m"
}
# restore.sh's output is kept in a temporary file to find the reason of a
# failure. out names that file only when mktemp created it; without one the
# output goes to /dev/null, which is never removed.
out=''
out=$(mktemp) || out=''
drop_out() { if [ -n "$out" ]; then rm -f -- "$out"; fi; }
# try_active: one restore of state/active; sets rc and reason.
try_active() {
    bash "$G/restore.sh" "${apply[@]}" 2>&1 | tee "${out:-/dev/null}"
    rc=${PIPESTATUS[0]}
    if [ -n "$out" ]; then
        reason=$(grep -m1 'FAIL' "$out" 2>/dev/null || tail -n 1 "$out" 2>/dev/null)
    else
        reason="restore.sh rc=$rc (no temporary file for its output: mktemp failed)"
    fi
}
active=$GRM_DEPLOY_ROOT/state/active
pair=unknown rc=1
if [ -f "$active" ] && [ ! -L "$active" ]; then
    pair="$(sed -n 's/^label=//p' "$active")/$(sed -n 's/^publication=//p' "$active")"
    say "scripts $(cat "${G%/scripts/deployment}/COMMIT" 2>/dev/null || echo unknown); state/active records $pair (commit $(sed -n 's/^commit=//p' "$active"), revision $(sed -n 's/^revision=//p' "$active"), $(sed -n 's/^datasets=//p' "$active") datasets)"
    try_active
    if [ "$rc" != 0 ]; then
        if ((dry)); then
            say "restore.sh for state/active failed (rc=$rc); dry run: would wait ${BOOT_RETRY_WAIT}s, then retry it once"
        else
            say "restore.sh for state/active failed (rc=$rc); retrying it once in ${BOOT_RETRY_WAIT}s"
            sleep "$BOOT_RETRY_WAIT"
        fi
        try_active
    fi
else
    reason="$active is missing or not a plain file"
fi
if [ "$rc" = 0 ]; then
    say "restore.sh ${apply[*]} for state/active passed"
    drop_out
    exit 0
fi
case $reason in 20[0-9][0-9]-*) reason=${reason#* } ;; esac
say "FAIL active=$pair reason=${reason:-restore.sh rc=$rc}"
mark "FAIL active=$pair rc=$rc reason=${reason:-unknown}"
if [ "$pair" = "$FB_LABEL/$FB_PUB" ]; then
    say "the failed pair is the fallback pair; nothing else to try"
    mark "FALLBACK not tried: the failed pair is $FB_LABEL/$FB_PUB"
    drop_out
    exit 1
fi
bash "$G/restore.sh" "$FB_LABEL" "$FB_PUB" "${apply[@]}"
frc=$?
say "FALLBACK to $FB_LABEL/$FB_PUB rc=$frc"
mark "FALLBACK to $FB_LABEL/$FB_PUB rc=$frc"
drop_out
[ "$frc" = 0 ] && exit 3
exit 1
