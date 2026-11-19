#!/usr/bin/env bash
# Run the pinned wrangler against the account with the token from
# GRM_SITE_TOKEN_FILE in its environment only, for the steps the release
# scripts do not cover (r2 bucket create, deploy, whoami). Logs the
# arguments, never the token. usage: GRM_SITE_TARGET=remote wrangler.sh ARGS...
source "$(dirname -- "${BASH_SOURCE[0]}")/site-lib.sh"
[ "$GRM_SITE_TARGET" = remote ] || die "wrangler.sh is for GRM_SITE_TARGET=remote" 2
APPLY=1
say "wrangler $*"
rc=0
(cd "$WORKER" && CLOUDFLARE_ACCOUNT_ID=$GRM_SITE_ACCOUNT_ID CLOUDFLARE_API_TOKEN=$(<"$GRM_SITE_TOKEN_FILE") "$WRANGLER" "$@") || rc=$?
say "wrangler $1 rc=$rc"
exit $rc
