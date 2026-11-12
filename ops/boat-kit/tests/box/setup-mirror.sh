#!/usr/bin/env bash
# KIT10 round 2 test box: lay out the pushed harness at the laptop paths the path-pinned tests read (PORT's
# drv-gates-v4.sh, KIT9's stub-regress, the drive-l2 test's setup), all on this box and all synthetic. PORT's
# synthetic master is made here by box/gen-port-master.sh from PORT's gen-root.py. Round 2 pushes no data set: no
# restore sets, no revision 12 curated/ or catalog/, no MTB source (the restore step is skipped this round).
set -eu
L=/home/layth/review H=/home/user/k10
[ -e /home/user/.mirror-done ] && { echo "mirror already set up"; exit 0; }
sudo -n mkdir -p "$L"; sudo -n chown -R "$(id -u):$(id -g)" /home/layth
P=$L/grm-port-r11-20261105
mkdir -p "$P/dryrun" "$P/drive" "$P/synth" "$L/grm-kit9-20261106/live"
cp -a "$H/port/dryrun/." "$P/dryrun/"; cp -a "$H/port/kit" "$P/kit"; cp -a "$H/live/drive.sh" "$P/drive/drive.sh"
cp -a "$H/live/drive.sh" "$L/grm-kit9-20261106/live/drive.sh"
bash "$H/box/gen-port-master.sh" "$H/port/kit/synth/gen-root.py" "$P/synth/master" "$P/synth/master.sha256"
sha256sum "$P/synth/master.sha256" "$L/grm-kit9-20261106/live/drive.sh" "$P/drive/drive.sh" "$P/dryrun/drv-gates-v4.sh" "$P/dryrun/drv-stub-runner-v4.sh"
(cd "$P/synth/master" && sha256sum -c --strict --quiet ../master.sha256) && echo "synth master ok"
touch /home/user/.mirror-done
