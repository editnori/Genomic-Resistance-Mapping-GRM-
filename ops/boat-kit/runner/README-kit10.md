# grm-boat.sh, KIT10 changes

These are KIT10's changes to grm-boat.sh (next/grm-boat.sh.next) next to the live runner (45fb0c9f…). The usage block at the top of the file has the full contract. tests/test-runner-kit10.sh checks each change below, and the check ids are in brackets.

## extend

`grm-boat.sh extend <id> --ttl SECONDS` sets the box's remaining lifetime to SECONDS from now. Boat's own extend sets that lifetime and doesn't add to it, so the runner refuses a TTL that would end the box earlier than its current archiveAfter. That refusal exits 2 and changes nothing [X1].

The runner reads archiveAfter before and after the call, and the exit is 86 in three cases [X5, X6, X7]:

- archiveAfter moved back
- archiveAfter didn't move
- archiveAfter ends more than 120 s short of now + SECONDS, with now taken before the call

A box that ends 60 s short passes [X8].

The last line on success is `END <epoch seconds>`, the box's new end [X2b]. drive.sh reads that line and sets its deadline from it. With no END line, the drive stops with 86 and closes the box (drive test H14c). When the box ends short of the ask but inside the 120 s, the drive's deadline is the END the box got [H14b].

## SigDrug vocabulary and the combine root

push-data takes GRM's SigDrug vocabulary copy only as `vocabulary/<name>.tsv` at the manifest's top level. It has to be all 13 of grm-r9-curate/vocabulary's files, each with its pinned sha256 [V0, V1, V2, V4]. A file named like one of the 13 is refused outside vocabulary/, whatever its folder or bytes, and the name is compared in lower case [W4, W5, W6]. That also closes tables/brand_generic_folds.tsv, which the live runner passed.

Two classes pass only when `--dest` is exactly /home/user/grm-combine, the drive's combine push:

- each part's pulled H-1 WORK, at `parts/<part>/h1work/...`, as text (json, jsonl, tsv, txt or md) [H1, H2b]
- one `MASTER-MANIFEST*.sha256` at the top [W7]

Any other destination refuses both with 97, and so does any other place in the tree: h1work/ outside parts/ [W1, W2, W10], a master manifest below the top [W3, W11], and the same combine root sent elsewhere [W8, W9]. Both classes are text receipts, at most 64 MiB each.

## Per-call staging

Each push, push-data and pull call stages its tar in its own folder, `stage-<lane>/call-<what>-XXXXXXXX`, and removes only that folder when it ends. Before KIT10, each call removed the whole `stage-<lane>`, so parallel calls in one lane deleted each other's tars [S1, S2].

The lane folder stays, even when it's empty [S3b]. Round 1 removed it when empty, and a call that had just made it then failed its mktemp with 86. S4 runs the runner's own stage_dir and unstage as two concurrent calls in one lane for 3000 rounds and needs zero failures.

## Not in this runner

KIT9's memory beacon and grm-memread.exe aren't in it, per main's ruling of 7 Nov [B1-B3, E1, E2 fail on KIT9's runner, as its control expects].

A file named MASTER-MANIFEST*.sha256 under audit/, receipts/ or tables/ still passes as an ordinary receipt for any destination, as on the live runner. It isn't a SigDrug path.
