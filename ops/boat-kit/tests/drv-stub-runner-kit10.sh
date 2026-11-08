#!/usr/bin/env bash
# KIT10's stub of grm-boat.sh for drive.sh's KIND=h1 (from REVIEW-PORT2's drv-stub-runner-v4.sh): logs each call to
# $SR_LOG with the GRM_SCP_FLOOR_MBPS it sees, keeps live boxes in $SR_LIVE, and answers like a box whose H-1 stage
# goes well. Knobs: SR_HOST fast|slow ; SR_SHA ; SR_STATUS done|running ; SR_UNION (a file of genome ids the box's
# union step reports) ; SR_PUSHDATA_FAIL=meta|vocab|fasta|combine (that push-data exits 97). KIT10 round 2:
# extend prints "END <now + ttl - SR_EXTEND_SHORT>" as grm-boat.sh extend does, or no END with SR_EXTEND_NOEND=1;
# SR_STATUS=paused answers as r11-box.sh does when an H-1 stage's guard holds it stopped; SR_RUNNING_POLLS=n answers
# running for n status polls, then SR_STATUS.
# KIT11: stage-local logs like push-data (its manifest as $SR_LOG.sl-<dest>, its --from as $SR_LOG.slfrom-<dest>);
# SR_STAGELOCAL_FAIL=meta|vocab|fasta makes that one exit 97.
L=${SR_LOG:?}; echo "$* | floor=${GRM_SCP_FLOOR_MBPS:-unset}" >> "$L"
[ "$1" = --notes ] && shift 2
c=$1; shift
now() { date -u +%FT%T.000Z; }
case "$c" in
  preflight) echo "stub preflight"; exit 0;;
  new)
    name=grm-$1; id=bx_stub$RANDOM$RANDOM; echo "$id $name $(now)" >> "$SR_LIVE"; echo "$id $name $*" >> "$SR_CREATED"
    mkdir -p "${GRM_BOAT_STATE:?}"; echo "$id $name $(date -u +%FT%TZ)" >> "$GRM_BOAT_STATE/created"; echo "$id";;
  push) exit 0;;
  push-data)
    d=""; m=""; while [ $# -gt 0 ]; do case "$1" in --dest) d=$2; shift;; --manifest) m=$2; shift;; esac; shift; done
    cp "$m" "$SR_LOG.pd-${d##*/}"
    case "${SR_PUSHDATA_FAIL:-}:$d" in meta:*grm-meta|vocab:*grm-vocab|fasta:*grm-fasta|combine:*grm-combine) exit 97;; esac; exit 0;;
  stage-local)
    d=""; m=""; f=""; while [ $# -gt 0 ]; do case "$1" in --dest) d=$2; shift;; --manifest) m=$2; shift;; --from) f=$2; shift;; esac; shift; done
    cp "$m" "$SR_LOG.sl-${d##*/}"; echo "$f" > "$SR_LOG.slfrom-${d##*/}"
    case "${SR_STAGELOCAL_FAIL:-}:$d" in meta:*grm-meta|vocab:*grm-vocab|fasta:*grm-fasta) exit 97;; esac; exit 0;;
  extend) ttl=""; while [ $# -gt 0 ]; do [ "$1" = --ttl ] && ttl=$2; shift; done
    echo "extended (stub)"; [ -n "${SR_EXTEND_NOEND:-}" ] || echo "END $(( $(date +%s) + ttl - ${SR_EXTEND_SHORT:-0} ))"; exit 0;;
  close) sed -i "/^$1 /d" "$SR_LIVE"; echo "closed $1" >> "$SR_CLOSED"; exit 0;;
  pull) # pull ID REMOTE LOCAL
    l=$(basename "$2"); d=$3/$l
    if [ "$l" = h1 ]; then mkdir -p "$d"; cp "${SR_UNION:?}" "$d/union.txt"; exit 0; fi
    mkdir -p "$d/run/receipts/h1" "$d/h1work/arms/default/ds1" "$d/h1work/j" "$d/h1work/k"
    echo "$(date -u +%FT%TZ) done" > "$d/run/status"
    for f in spec.json inputs.json summary.json species-groups.json raw-manifest.txt j/contrasts.json k/manifest.json arms/default/ds1/runs.json; do echo '{}' > "$d/h1work/$f"; done
    echo '{"work_blocks_bytes": 1, "all_roots_blocks_bytes": 1, "arms_blocks_bytes": 1, "pool": null, "peak_rss_bytes": 1}' > "$d/run/receipts/h1/disk-measured.json"
    (cd "$d" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS); exit 0;;
  run)
    log=""; while [ $# -gt 0 ]; do case "$1" in --log) log=$2; shift 2;; --cwd) shift 2;; --) shift; break;; bx_*) shift;; *) break;; esac; done
    case "$log" in
      stage-env|fasta-pin) printf '%s\n' "$*" >> "$SR_LOG.env";;
      host) echo "cpu=Stub CPU"; echo "${SR_SHA:-$(printf 'a%.0s' {1..64})}  /home/user/grm/target/release/grm"; echo "HOST ${SR_HOST:-fast}";;
      prep) echo "numpy 1.26.4";;
      union) echo "UNION $(wc -l < "$SR_UNION") $(sha256sum "$SR_UNION" | cut -d' ' -f1)";;
      verify) echo "RAW 9 files"; echo "revision-12 gate: PASSED (stub)";;
      status)
        n=$(( $(cat "$L.polls" 2> /dev/null || echo 0) + 1 )); echo $n > "$L.polls"
        s=${SR_STATUS:-done}; [ -n "${SR_RUNNING_POLLS:-}" ] && [ "$n" -le "$SR_RUNNING_POLLS" ] && s=running
        if [ "$s" = paused ]; then echo "STATE paused | guard disk: 30064771072 B free, below DISK_PAUSE_GIB 30 GiB (paused 2026-11-07T12:00:00Z disk_free=30064771072 eff_mem=99999999999) | running stage=stub"
        else echo "STATE $s | stub"; fi;;
      ckpt-*) echo "CKPT stub";;
    esac; exit 0;;
  *) echo "stub runner: $c?" >&2; exit 2;;
esac
