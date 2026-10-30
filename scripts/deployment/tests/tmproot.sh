# Sourced by the deployment test suites, which delete their temporary root on
# exit. make_test_root prints a fresh, empty temporary directory or fails.
# The obvious T=$(cd -- "$(mktemp -d)" && pwd -P) is unsafe: when mktemp fails
# (TMPDIR missing, say) it prints nothing, `cd -- ""` succeeds, and the suite
# takes the current directory as its root and deletes it on exit. Run from the
# repository, that deleted the whole working tree.
make_test_root() {
    local t repo
    t=$(mktemp -d) || { echo "cannot create a temporary directory (TMPDIR=${TMPDIR-unset})" >&2; return 1; }
    t=$(cd -- "$t" && pwd -P) || return 1
    repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd -P)
    case $t/ in
        /*/*/) ;;
        *) echo "refusing temporary directory '$t': too close to /" >&2; return 1 ;;
    esac
    case $repo/ in
        "$t"/*) echo "refusing temporary directory '$t': it contains the repository" >&2; return 1 ;;
    esac
    [ -z "$(ls -A -- "$t")" ] || { echo "refusing temporary directory '$t': not empty" >&2; return 1; }
    printf '%s\n' "$t"
}
