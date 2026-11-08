# KIT10 install

install-kit10.sh installs KIT9's and KIT10's changes on the laptop in one run. KIT10's drive.sh.next and grm-boat.sh.next already hold KIT9's changes, with REVIEW-KIT9's fixes and without KIT9's beacon and grm-memread.exe (main's ruling of 7 Nov), so KIT9's own drive.sh.next and grm-boat.sh.next don't get installed. This is round 2's set, with REVIEW-KIT10's fixes. It replaces round 1's set, which was never installed.

## Files

Each new file replaces the live copy that KIT10's tests ran against, which is k10/live for the drive and runner and PORT's kit for the kit files. Three kit files are new. The mode comes from the file it replaces, and a new kit file takes the mode of r11-box.sh or stage-manifest.py.

    next/drive.sh.next               -> /home/layth/review/grm-port-r11-20261105/drive/drive.sh
      live  db62192168156fc773c46fe101c188e3cfa2192312cb46d356a72a81f2d860f7
      new   14ef9a19fa7af2f0a5e0ce00998a2b8b1d377f1ccce4511ab20deb47a77e27d2
    next/grm-boat.sh.next            -> /home/layth/review/grm-boat-runner-20261031/grm-boat.sh
      live  45fb0c9ff27c7a851f13c5bcbfaa1c99388a47d61f7511ce409ba69ae05c3cd2
      new   dd3f4842c07fdf9ccc3a5da4f840066e3a1517aaea1403e4fe3ed752e9b2db5e
    next/kit/mtb-box.sh.next         -> /home/layth/review/grm-port-r11-20261105/kit/mtb-box.sh
      live  0a42c769c633435f7e1200522bfd2aaec53b61843935d225d56f0a8b4f850361
      new   f67db557c7e642efb01e4a01184cc1103424e97a16e787d2d4a997569ee5f072
    next/kit/r11-box.sh.next         -> /home/layth/review/grm-port-r11-20261105/kit/r11-box.sh
      live  4920a98e6066fe5b39778c293ebb43f74a7a74244065e1be6795e474100f8d56
      new   b618def77fce3d39f3d80f8ea9a6fe359850f5c77d4cdec330ff5e304b85bef1
    next/kit/stage-manifest.py.next  -> /home/layth/review/grm-port-r11-20261105/kit/stage-manifest.py
      live  36a9be604d0cc23c515f498a46ff5d1492ea5c6b8f01774cd842d542daf99699
      new   68cf7b4cd25fa4ddea2df407ddf51b118dda2a3dba6d2c8f5180d62c5e1eddc5
    next/kit/chain-lib.sh.next       -> /home/layth/review/grm-port-r11-20261105/kit/chain-lib.sh
      new   1e9514d0d3fdf705d043b3edb037bc889355871a76d2c96694601f1aafc28e3b
    next/kit/h1-box.sh.next          -> /home/layth/review/grm-port-r11-20261105/kit/h1-box.sh
      new   b6f32caf66bef5e12403e242ed5e2666da78828f4639212919be4b8afa98e268
    next/kit/h1-manifest.py.next     -> /home/layth/review/grm-port-r11-20261105/kit/h1-manifest.py
      new   12f93492cb482432960fd818e4cd843cf5d5f9ef86c4bba93326a1c6abe46e09

The three target roots come from the harness: the drive from drv-gates-v4.sh's default, the runner from restart-chain.sh, and the kit from make-harness.sh. If one is wrong, the hash check stops the install before anything changes. diffs/next-vs-live has each change as a unified diff, and diffs/next-vs-r1 has round 2's changes to round 1's set.

mtb-box.sh.next is byte-equal to round 1's (f67db557…). Its restore test isn't rerun this round, so its evidence is run 5's restore receipts (REPORT.md).

next/run-h1.conf, next/stages-h1.tsv and next/README-runner-kit10.md aren't installed. The first two are drafts of H-1's run file and stage table, and drive.sh refuses them while they hold placeholders. The README is the runner's notes for KIT10.

## Steps

install-kit10.sh ships beside this file and is covered by SHA256SUMS. Read it before the run.

Copy this handback dir to the laptop, cd into it, and run the script there when no drive.sh, grm-boat.sh or chain is running:

    bash install-kit10.sh

Step 1 checks everything before any change:

- the handback against SHA256SUMS
- that no bash is running drive.sh, grm-boat.sh or chain-tests.sh (an editor or a log viewer on one of them doesn't count)
- each new file's hash and syntax
- each live file's hash
- that each target's dir and this dir are writable
- that no pre-kit10/ dir is there yet
- that the kit repo doesn't gitignore a new kit file

If any check fails, it prints STOP and changes nothing.

Step 2 puts each file in through a temp file beside its target. The temp file takes the live file's mode with `chmod --reference` and replaces the target with `mv`, so a script already running keeps reading its old file.

Each file it replaces is kept in pre-kit10/ in this handback dir. That dir sits outside the kit repo on purpose. The runner's push sends `git ls-files -co --exclude-standard`, untracked files included, so a backup inside the kit would go to every box.

If anything in step 2 fails, the script rolls back what it changed and prints why. It removes the temp file, puts each replaced file back from pre-kit10/ through a temp file and `mv`, and removes each new file it made. It then checks every live hash. When all of them match, it removes pre-kit10/ and says "ROLLED BACK", so a rerun starts clean. If one doesn't match, it says "ROLLBACK INCOMPLETE" and keeps pre-kit10/ for the rollback below.

The last line is "KIT10 installed: 8 files, syntax ok". The kit repo then shows three changed files (mtb-box.sh, r11-box.sh and stage-manifest.py) and three untracked ones (chain-lib.sh, h1-box.sh and h1-manifest.py). A box push sends all six.

To put the live files back after an install, run this from the handback dir. Unlike the install, it copies over each live file in place, so run it only when no drive.sh, grm-boat.sh or chain is running. It checks the live hashes and removes pre-kit10/ only when they all match, so the install can run again:

    D=/home/layth/review/grm-port-r11-20261105 R=/home/layth/review/grm-boat-runner-20261031 P=$PWD/pre-kit10
    cp -p $P/drive.sh $D/drive/drive.sh; cp -p $P/grm-boat.sh $R/grm-boat.sh
    for f in mtb-box.sh r11-box.sh stage-manifest.py; do cp -p $P/$f $D/kit/$f; done
    rm -f $D/kit/chain-lib.sh $D/kit/h1-box.sh $D/kit/h1-manifest.py
    printf '%s\n' "db62192168156fc773c46fe101c188e3cfa2192312cb46d356a72a81f2d860f7  $D/drive/drive.sh" \
      "45fb0c9ff27c7a851f13c5bcbfaa1c99388a47d61f7511ce409ba69ae05c3cd2  $R/grm-boat.sh" \
      "0a42c769c633435f7e1200522bfd2aaec53b61843935d225d56f0a8b4f850361  $D/kit/mtb-box.sh" \
      "4920a98e6066fe5b39778c293ebb43f74a7a74244065e1be6795e474100f8d56  $D/kit/r11-box.sh" \
      "36a9be604d0cc23c515f498a46ff5d1492ea5c6b8f01774cd842d542daf99699  $D/kit/stage-manifest.py" | sha256sum -c --strict && rm -rf $P

tests/test-install-kit10.sh runs the script and this rollback on a fake tree with the live hashes and modes. logs/run6 has its run on Boat. It covers:

- a clean install
- this rollback, and an install after it
- a second run
- a drifted live file
- a changed handback file
- a read-only target dir
- a failure after the first file is in
- a failure after mktemp
- a failed move of the last file
- the running-script check

tests/test-restore-v7.sh.next and tests/test-side-re.sh.next replace main's copies of those tests. Those copies aren't in the harness, so this install leaves them out.
