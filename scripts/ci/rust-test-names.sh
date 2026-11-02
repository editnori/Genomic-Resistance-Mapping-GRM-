#!/usr/bin/env bash
# The names of the Rust tests that `cargo test --release --locked` builds, one
# "TARGET NAME" line per test, sorted. TARGET is the test's crate root
# (src/lib.rs, src/main.rs, tests/NAME.rs), or "doc" for a doc test, whose
# name drops its "(line N)" so that moving the doc test keeps its name. The
# list holds ignored tests too, each ending in " [ignored]". It must hold as
# many tests as libtest itself counts, or the script fails.
#
#   rust-test-names.sh           print the list
#   rust-test-names.sh --check   fail unless it equals tests/rust-test-names.txt
#   rust-test-names.sh --write   rewrite tests/rust-test-names.txt
#
# docs/ci.md says when to rewrite the file.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.."
FILE=tests/rust-test-names.txt
HEADER='# Rust test names: scripts/ci/rust-test-names.sh --write regenerates this file (docs/ci.md).'

# Prints "TARGET NAME" for each test that `cargo test -- --list "$@"` lists.
# Cargo names a target on stderr before that target lists its tests on stdout.
list() {
  CARGO_TERM_COLOR=never cargo test --release --locked -- --list "$@" 2>&1 | awk '
    $1 == "Running" { t = ($2 == "unittests") ? $3 : $2; next }
    $1 == "Doc-tests" { t = "doc"; next }
    /: test$/ {
      if (t == "") { print "rust-test-names.sh: a test before any target: " $0 > "/dev/stderr"; exit 3 }
      sub(/: test$/, ""); print t " " $0; n++; next
    }
    /^[0-9]+ tests?, [0-9]+ benchmarks?$/ { counted += $1 }
    END {
      if (n != counted) { print "rust-test-names.sh: read " n " tests, but libtest counted " counted > "/dev/stderr"; exit 3 }
    }'
}

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
list > "$tmp/all"
list --ignored > "$tmp/ignored"
# Mark the ignored tests, then drop the doc tests' line numbers.
awk 'FILENAME == ARGV[1] { ignored[$0]; next } { print $0 (($0 in ignored) ? " [ignored]" : "") }' "$tmp/ignored" "$tmp/all" |
  sed -E '/^doc /s/ \(line [0-9]+\)( \[ignored\])?$/\1/' | LC_ALL=C sort > "$tmp/names"
if [ "$(grep -c ' \[ignored\]$' "$tmp/names")" != "$(wc -l < "$tmp/ignored")" ]; then
  echo "rust-test-names.sh: an ignored test is missing from the full list" >&2
  exit 3
fi
{ echo "$HEADER"; cat "$tmp/names"; } > "$tmp/file"
summary="$(wc -l < "$tmp/names") tests, $(grep -c ' \[ignored\]$' "$tmp/names" || true) of them ignored"

case ${1-} in
  "") cat "$tmp/names" ;;
  --write)
    cp "$tmp/file" "$FILE"
    echo "wrote $FILE: $summary"
    ;;
  --check)
    if diff -u --label "$FILE" --label "cargo test -- --list" "$FILE" "$tmp/file"; then
      echo "$FILE matches: $summary"
    else
      echo "The Rust test names differ from $FILE: lines with - left the build, lines with + are new ($summary)." >&2
      echo "If the change is meant, run scripts/ci/rust-test-names.sh --write and commit $FILE with it (docs/ci.md)." >&2
      exit 1
    fi
    ;;
  *)
    echo "usage: $0 [--check | --write]" >&2
    exit 2
    ;;
esac
