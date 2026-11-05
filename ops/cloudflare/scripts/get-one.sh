#!/usr/bin/env bash
# One object read for verify_objects (site-lib.sh): get-one.sh KEY FILE.
# rc 0 when written, 1 when the key does not exist.
source "$(dirname -- "${BASH_SOURCE[0]}")/site-lib.sh"
get "$1" "$2"
