#!/usr/bin/env bash
# KIT10 round 2 (controller): build the test box's harness as one tarball from round 2's work tree and handback, plus
# round 1's set for the round 1 controls, all read from this controller. No data set goes in; PORT's synthetic master
# is made on the box by box/gen-port-master.sh. Two round 1 scripts stay out: gen-synth.sh, the H-1 dry run's
# generator (the dry run is out of scope this round), and restart-chain.sh, which drives the box from the laptop.
# Round 3 (REVIEW-KIT10-R2 L1): the handback is round 3's, from its own work tree, and round 2's handback goes in as
# r2/handback for the install's round 2 control.
# usage: make-harness.sh OUT_DIR   (writes OUT_DIR/k10.tar.gz, its .sha256 and k10-files.sha256, every file's hash)
set -eu
R3=/home/user/grm-kit10-r3 R2=/home/user/grm-kit10-r2 R1=/home/user/grm-kit10 K=/home/user/grm-kit10-r3/work/k10
O=${1:?out dir} B=${1:?out dir}/k10
rm -rf "$B" "$O/k10.tar.gz" "$O/k10.tar.gz.sha256"; mkdir -p "$B/box" "$B/handback" "$B/r1/handback" "$B/r2/handback"
cp -a "$K/cand" "$K/kit" "$K/kit9" "$K/live" "$K/port" "$K/prov" "$K/tests" "$K/drafts" "$B/"
for f in "$K"/box/*; do case "${f##*/}" in gen-synth.sh|restart-chain.sh) ;; *) cp -a "$f" "$B/box/";; esac; done
cp -a "$R3/install-kit10.sh" "$R3/INSTALL.md" "$R3/next" "$B/handback/"
cp -a "$R2/install-kit10.sh" "$R2/INSTALL.md" "$R2/next" "$B/r2/handback/"
cp -a "$R1/work/k10/cand" "$R1/work/k10/kit" "$B/r1/"
cp -a "$R1/install-kit10.sh" "$R1/INSTALL.md" "$R1/next" "$B/r1/handback/"
chmod +x "$B"/tests/*.sh "$B"/box/*.sh "$B"/cand/*.sh "$B"/kit/*.sh
(cd "$O" && find k10 -type f | LC_ALL=C sort | xargs sha256sum) > "$O/k10-files.sha256"
tar -C "$O" --sort=name -czf "$O/k10.tar.gz" k10
(cd "$O" && sha256sum k10.tar.gz > k10.tar.gz.sha256)
echo "harness $(wc -l < "$O/k10-files.sha256") files, $(stat -c %s "$O/k10.tar.gz") B, $(cat "$O/k10.tar.gz.sha256")"
