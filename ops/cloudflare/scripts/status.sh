#!/usr/bin/env bash
# Print state/active and state/previous, and the last lines of the log.
# Read-only. usage: status.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/site-lib.sh"
f=$(mktemp)
for which in active previous; do
    if get "state/$which" "$f"; then printf "%-8s %s" "$which" "$(cat "$f")"; echo; else echo "$which (none)"; fi
done
rm -f "$f"
if [ -f "$GRM_SITE_LOG" ]; then echo "log $GRM_SITE_LOG:"; tail -n 5 "$GRM_SITE_LOG"; fi
