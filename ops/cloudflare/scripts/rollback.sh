#!/usr/bin/env bash
# Switch the site back to the pair state/previous names, with every check
# activate.sh makes. A second rollback switches forward again.
# usage: rollback.sh [--apply]
source "$(dirname -- "${BASH_SOURCE[0]}")/site-lib.sh"
parse_args "$@"
[ ${#ARGS[@]} -eq 0 ] || die "usage: rollback.sh [--apply]" 2
f=$(mktemp)
get state/previous "$f" || { rm -f "$f"; die "state/previous is missing: nothing to roll back to" 2; }
set -- $(read_pair "$f") || die "state/previous is not a valid pair record"
rm -f "$f"
say "rollback target from state/previous: $1 with $2"
if ((APPLY)); then exec bash "$HERE/activate.sh" --rollback "$1" "$2" --apply; fi
exec bash "$HERE/activate.sh" --rollback "$1" "$2"
