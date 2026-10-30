#!/usr/bin/env bash
# Boot hook (run from /root/onstart.sh): bring back the pair state/active
# records on the live port, through the deployment toolkit. restore.sh with no
# pair restores state/active: it verifies the release and the publication
# against their manifests, the approval and state/active itself, resolves an
# interrupted start, starts the pair as a supervised start, requires it to own
# the port and pass smoke.sh, and records it. When that pair already serves and
# is recorded, it restarts nothing.
#
# Fallback (main's ruling, 2026-10-29; moved to revision 8, 2026-10-30): when
# that fails closed, restore.sh ui-r8-690815a rev8-273-20261029 --apply, through
# the same checks, never
# serve-live.sh. Both outcomes are logged (FAIL active=..., FALLBACK ... rc=...)
# and recorded in BOOT-FALLBACK in the root. If the fallback fails too, nothing
# serves.
# usage: boot-active.sh [--dry-run]   (dry run: restore.sh without --apply,
#   and BOOT-FALLBACK is not written)
# BOOT_DEPLOY_ROOT and BOOT_SCRIPTS replace the live root and scripts (tests).
# Exit 0 when the active pair serves, 3 when the fallback serves, else 1.
set -u
export GRM_DEPLOY_ROOT=${BOOT_DEPLOY_ROOT:-/root/grm-deploy-20261025}
G=${BOOT_SCRIPTS:-/root/grm-ops-20261029/scripts/deployment}
FB_LABEL=ui-r8-690815a FB_PUB=rev8-273-20261029
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
out=$(mktemp) || out=/dev/null
active=$GRM_DEPLOY_ROOT/state/active
pair=unknown rc=1
if [ -f "$active" ] && [ ! -L "$active" ]; then
    pair="$(sed -n 's/^label=//p' "$active")/$(sed -n 's/^publication=//p' "$active")"
    say "scripts $(cat "${G%/scripts/deployment}/COMMIT" 2>/dev/null || echo unknown); state/active records $pair (commit $(sed -n 's/^commit=//p' "$active"), revision $(sed -n 's/^revision=//p' "$active"), $(sed -n 's/^datasets=//p' "$active") datasets)"
    bash "$G/restore.sh" "${apply[@]}" 2>&1 | tee "$out"
    rc=${PIPESTATUS[0]}
    reason=$(grep -m1 'FAIL' "$out" 2>/dev/null || tail -n 1 "$out" 2>/dev/null)
else
    reason="$active is missing or not a plain file"
fi
if [ "$rc" = 0 ]; then
    say "restore.sh ${apply[*]} for state/active passed"
    rm -f -- "$out"
    exit 0
fi
case $reason in 20[0-9][0-9]-*) reason=${reason#* } ;; esac
say "FAIL active=$pair reason=${reason:-restore.sh rc=$rc}"
mark "FAIL active=$pair rc=$rc reason=${reason:-unknown}"
if [ "$pair" = "$FB_LABEL/$FB_PUB" ]; then
    say "the failed pair is the fallback pair; nothing else to try"
    mark "FALLBACK not tried: the failed pair is $FB_LABEL/$FB_PUB"
    rm -f -- "$out"
    exit 1
fi
bash "$G/restore.sh" "$FB_LABEL" "$FB_PUB" "${apply[@]}"
frc=$?
say "FALLBACK to $FB_LABEL/$FB_PUB rc=$frc"
mark "FALLBACK to $FB_LABEL/$FB_PUB rc=$frc"
rm -f -- "$out"
[ "$frc" = 0 ] && exit 3
exit 1
