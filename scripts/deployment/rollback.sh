#!/usr/bin/env bash
# Switch the live service back to the release and publication that
# activate.sh recorded in state/previous, with the same checks, approval and
# automatic restore as any activation. A second rollback switches forward again.
# usage: rollback.sh [--apply]
# The target must verify; the live pair need not. When the live pair fails its
# checks (a damaged release, say), it is replaced only if it is a managed grm
# serve of exactly the pair state/active records, and it is not restarted if
# the target fails (activate.sh --rollback). It needs a running managed
# service; when nothing serves the port, start the wanted pair with
# restore.sh LABEL PUBLICATION instead.
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/common.sh"
parse_args "$@"
[ ${#ARGS[@]} -eq 0 ] || die "usage: rollback.sh [--apply]" 2
need python3 sha256sum
verify_state previous
say "rollback target from state/previous: $ST_LABEL with $ST_PUB"
if ((APPLY)); then
    exec bash "$HERE/activate.sh" --rollback "$ST_LABEL" "$ST_PUB" --apply
fi
exec bash "$HERE/activate.sh" --rollback "$ST_LABEL" "$ST_PUB"
