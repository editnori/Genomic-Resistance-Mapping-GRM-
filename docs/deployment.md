# Deployment

genomicresistance.com is served by one `grm serve` process on a single Vast container, bound to 127.0.0.1:8641, and published through a Cloudflare Tunnel whose connector runs on the same container. The scripts in `scripts/deployment/` stage, switch, roll back and restore that process from immutable builds and report copies.

They never deploy by default. `publish.sh`, `build-release.sh`, `stage.sh`, `activate.sh`, `rollback.sh`, `restore.sh` and `approve-legacy.sh` change the deployment root or the service only with `--apply`; without it they run their read-only checks, print what they would do, and write nothing under the root. `connector.sh` has no dry run: `check` changes nothing, but `run` and `stop` act at once (see [Tunnel connector](#tunnel-connector)). `smoke.sh` and `status.sh` only send GET requests to 127.0.0.1, and they write nothing but a temporary directory that `smoke.sh` creates with `mktemp` and removes.

## Limitations

- One host, one origin process. Activation stops the old process before starting the new one (the port allows one listener), so the site returns errors for the second or two between them. There is no standby host; if the container or its disk is lost, rebuild it from a source commit and a verified report backup.
- The container has no service manager. Nothing restarts `grm serve` when it exits; `restore.sh` must be run after a crash or a container restart. `/root/onstart.sh` restarts only the tunnel connector at boot; see [After a container restart](#after-a-container-restart).
- The tunnel, its public hostname route, the DNS record and the Cloudflare Access bypass application live in the Cloudflare account and are not managed here. The scripts check only the local origin. The public URL depends on Cloudflare, the domain registration and the connector's outbound HTTP/2 connections (QUIC failed from the Vast host).
- The scripts were tested by a per-file `bash -n`, by `tests/dry-run.sh`, and by `tests/apply.sh`, which runs every `--apply` path in a temporary root against a stub server compiled from C (not `grm`) and a stub cloudflared (a shell script that execs `sleep`) (see [Checks](#checks)), and by a mutation harness that reintroduces each reviewed defect (`data/work/deployment-mutation-20261027/`). None of them has been run against the live service, the real `grm` or the real cloudflared, and none checks the public site. Signals in the tests come from test hooks at named steps of the scripts and of the forked child before it execs, plus external TERM and KILL; a real `grm` or cloudflared that takes longer than about 10 s to listen or exit has not been tested.
- A started process is tracked by a lock on a spawn marker that it inherits at the fork (see [Supervised starts](#supervised-starts)). That relies on the server keeping inherited descriptors open, which the stub does and Rust's standard library and Go's runtime do (they mark only their own descriptors close-on-exec); it has not been checked against the real binaries. A server that closed the marker would still be found by its registered PID and start time and its process group, but a process it forks that then closes the marker and leaves the group would not be.
- `stage.sh` takes no root lock. A second stage on the same spare port is refused while the first one's journal names a live owner, but two stages started at the same instant on one port are not supported.
- A process that survives SIGKILL (stuck in the kernel) cannot be stopped: activation then exits 4 without restarting the old pair, and `restore.sh` refuses, both keeping `serve-PORT.starting` and the marker so the start stays identified. That path is not exercised by the tests.
- A survivor the scripts cannot see (another user's process, or one that made itself non-dumpable, so its descriptors are unreadable) still holds the marker's lock: `restore.sh` and `stage.sh` then keep the journal, start nothing on that port and fail (tested), and the operator must find it, for example as root with `fuser` on the marker named in the journal. The same holds when the marker file was deleted: the scripts cannot prove that nothing survives, so they keep the journal (tested); remove `serve-PORT.starting` and `spawn/NAME` by hand only after checking.
- A rollback restart that cannot write its start journal or marker (a full disk) still starts the old pair; without a journal its child does not wait for a registration, and without a marker it is found only as the script's own child, so if the script is then killed outright, that server is found only as the live listener. Not tested.
- Files are opened only after a check that they are absent or regular and not links, and every descriptor is checked again after the open, so a link, FIFO or special file swapped in between is refused without blocking. Two gaps remain: a dangling link swapped in at that moment can have its target created empty before the check on the descriptor refuses it, and tree verification (`sha256sum -c` over a release or report tree, and `publish.sh` reading a source manifest) checks for special files before reading, not on the descriptor. Both need someone racing the scripts with write access to the root, which could change anything anyway.
- Approvals, legacy `GO-LABEL` files and state must be strict printable-ASCII records (see [Approvals and state](#approvals-and-state)). A hand-written `GO-LABEL` with tabs, CRLF endings, extra lines or other characters in its values is refused; stage that pair instead.

## Configuration and layout

Configuration comes from the environment, never from the scripts:

| Variable | Meaning |
|---|---|
| `GRM_DEPLOY_ROOT` | Absolute deployment root (required) |
| `GRM_PORT` | Live loopback port, default 8641 |
| `GRM_EXPECT_SCHEMA`, `GRM_EXPECT_REVISION`, `GRM_EXPECT_DATASETS` | What `publish.sh`, `stage.sh` and `approve-legacy.sh` require of a report index: default 2, 7 and 237 |

The approval that `stage.sh` or `approve-legacy.sh` writes records the schema, revision and dataset count it checked. `activate.sh`, `rollback.sh`, `restore.sh` and `status.sh` check each pair against its own approval and ignore `GRM_EXPECT_*` (they say so when it is set), so a 275-dataset candidate and a 237-dataset rollback target need no shared setting.

Under the root:

| Path | Contents |
|---|---|
| `published/NAME/` | `report/` and `SOURCE.sha256`, read-only |
| `release/LABEL/` | `grm`, `dist/`, `COMMIT`, `MANIFEST.sha256`, read-only |
| `src-LABEL/`, `src-LABEL.commit` | The archived source a release was built from |
| `stage/LABEL--NAME--TIME.XXXXXX/` | Candidate check receipts (unique per run) |
| `go/LABEL--NAME` | Approval of one build of a release with one copy of a publication |
| `state/active`, `state/previous` | The pair that should be serving, and the rollback target |
| `serve-PORT.starting` | Journal of a server start on PORT that has not finished; see [Supervised starts](#supervised-starts) |
| `serve-PORT.pid` | `pid=N start=N label=L publication=P [token=T mark=M mark_id=D:I]` of the committed server on PORT, read from the process itself; it identifies a damaged server for rollback |
| `spawn/TOKEN.XXXXXX` | Marker of one supervised start, held open (and locked) by the processes it started; see [Supervised starts](#supervised-starts) |
| `deploy.log`, `logs/serve-PORT.log`, `deploy.lock` | Operator log, server output and the root lock |

Labels and names may contain only letters, digits, `.`, `_` and single `-`.

A tree verification, used for every publication and release below, requires that the tree holds only directories and regular files (a symlink, device, FIFO or socket anywhere fails), that its manifest lists each file once with a matching SHA-256, and that nothing is missing or unlisted. Only a release's own `MANIFEST.sha256` is left out of its listing; a stray `MANIFEST.sha256` inside a report is an unlisted file. A release's manifest must list `grm`, `COMMIT` and `dist/index.html`, and `COMMIT` must be one full commit id.

Which script checks what:

| Script | Checks |
|---|---|
| `publish.sh` | The source tree against its manifest and its index against `GRM_EXPECT_*`, then the copy again |
| `stage.sh` | Release and publication trees, the index against `GRM_EXPECT_*`, then the served bytes |
| `activate.sh`, `restore.sh` | Release and publication trees, the approval, the index against the approval, then the served bytes. `activate.sh` checks the live pair the same way, approval included, before it stops anything |
| `rollback.sh` | `state/previous` against the current build, report copy and approval, then everything `activate.sh` checks, except that a live pair failing its checks is replaced when it is the managed pair `state/active` records (see [Roll back](#roll-back)) |
| `approve-legacy.sh` | Release and publication trees, the index against `GRM_EXPECT_*`, the legacy file, and the served bytes of the live pair |
| `status.sh` | No manifests: it reports the state files, `serve-PORT.pid`, any start journal, what listens, whether `state/active` records it, and runs `smoke.sh` against the live pair's approval (or against `GRM_EXPECT_*` when the live pair has no valid approval) |

`smoke.sh` compares the served `/`, `corpus.json`, the first dataset report, the app's `index-*.js` and `index.json` byte for byte with the release and publication, and checks the served index's schema, revision and dataset count.

### Approvals and state

`go/LABEL--NAME` is one line of unique `key=value` fields, each compared exactly: `stage=0`, `browser=pass|skipped`, `commit` (the release's `COMMIT`), `publication`, `release_manifest` and `report_manifest` (the SHA-256 of the release's `MANIFEST.sha256` and of the publication's `SOURCE.sha256`), `schema`, `revision`, `datasets`, and `receipts`. A missing, repeated, malformed or unknown field fails, so does a rebuilt release or a re-published report under the same name.

`state/active` and `state/previous` hold one field per line: the label, publication, commit, both manifest digests and the approval's schema, revision and dataset count, then the time. `rollback.sh` and an argument-free `restore.sh` refuse a state file that no longer matches the current build, report copy and approval.

Every record (approvals, legacy `GO-LABEL` files, state, `serve-PORT.starting`, `serve-PORT.pid`) is read whole, once, before any field is parsed. It must be a regular file (a link, directory or FIFO is refused without being opened) of at most 64 KiB, hold only printable ASCII and newlines, and end with a newline. A NUL, CR (so a CRLF line ending), tab, other control byte or non-ASCII byte anywhere refuses it, so nothing after a NUL can hide a conflicting field. Approvals and legacy files are exactly one line of fields separated by single spaces; state files are one field per line with no blank line. Values are limited to letters, digits and `._:/+-`.

Every record is written the same way: the target must be absent or a regular file (a link, directory, FIFO or other special file is refused before anything is written), the content goes to a new temporary file beside it, is compared byte for byte, is renamed over the target with `mv -T` (which fails on a directory rather than moving into it), and is compared byte for byte again. A state file is then parsed back and must hold exactly the expected fields. `activate.sh` writes `state/previous` and `state/active` this way before it stops anything, so a write that fails, or a state file that differs after its rename, refuses the activation with the old service untouched (exit 2). An interrupted write leaves only a hidden temporary file beside the old record, which nothing reads.

All file access goes through `scripts/deployment/lib.sh`, which `common.sh`, `connector.sh` and `smoke.sh` share. Records, digests and single lines (a release's `COMMIT`, `cloudflared.sha256`, `supervisor.pid`) are read by Python with `O_NOFOLLOW | O_NONBLOCK`, checking the descriptor it got. Descriptors the shell holds (`deploy.lock`, `connector.lock`, spawn markers, `deploy.log`, `logs/serve-PORT.log`) are opened only when the path is absent or a regular file, read-write (so a FIFO cannot block the open), and then checked on the descriptor itself: it must be a regular file at that path. So a FIFO `deploy.log` is skipped, a special `deploy.lock`, server log, `spawn/` or marker refuses the start, and nothing blocks.

### The live root

The live root on Vast is `/root/grm-deploy-20261025` with publication `rev7-20261025`. Releases there were approved before this toolkit by hand-written `GO-LABEL` files, which name a commit but no publication or report copy. The scripts no longer accept them. `approve-legacy.sh LABEL NAME` migrates one for the pair that serves the port now, after checking the release and publication trees, the index against `GRM_EXPECT_*`, that `GO-LABEL` is one line of unique fields that passes the byte checks in [Approvals and state](#approvals-and-state) and says `stage=0`, names the release's commit and names no other publication, and the live health check. It writes `go/LABEL--NAME` with `browser=legacy` and the legacy file's SHA-256. Any other pair, such as an older release kept as a rollback target, needs `stage.sh`. That root has no `state/` yet, so `rollback.sh` refuses until the first `activate.sh`, and the first `restore.sh` there needs an explicit pair.

Because activation needs an approval for the live pair (it is the rollback target), migrate it first:

```bash
bash $G/approve-legacy.sh ui-r2-9f193bf rev7-20261025             # dry run
bash $G/approve-legacy.sh ui-r2-9f193bf rev7-20261025 --apply
bash $G/restore.sh ui-r2-9f193bf rev7-20261025 --apply            # records state/active; restarts nothing
```

Copy the scripts from a reviewed commit to the host; they find each other by their own directory, so any location works:

```bash
git archive HEAD scripts/deployment | ssh vast 'mkdir -p /root/grm-ops && tar -x -C /root/grm-ops'
# on the host
export GRM_DEPLOY_ROOT=/root/grm-deploy-20261025
G=/root/grm-ops/scripts/deployment
bash $G/status.sh                   # read-only: state, live release, live health
```

## Local preview

```bash
cargo build --release --locked
npm --prefix web ci && npm --prefix web run build
./target/release/grm --data data/backups/revision7-published-20261025 serve --port 8641
```

`--data` names a root that contains `report/`; the verified backup has that layout. Open <http://127.0.0.1:8641>. For the Vite dev server, set `GRM_REPORT_DIR=data/backups/revision7-published-20261025/report`.

## Verified backup

`data/backups/revision7-published-20261025` holds the published revision-7 report (252 files: index, corpus, 237 dataset and 13 species reports) with `SOURCE.sha256`, plus the curated inputs, logs and receipts that validate it; its `README.md` describes the checks. Do not modify it. Verify it before use:

```bash
B=data/backups/revision7-published-20261025
(cd $B/report && sha256sum --quiet --strict -c ../SOURCE.sha256)
python3 -I -B $B/validation/run-check-results.py
```

To make a publication from it on a host, copy `report/` and `SOURCE.sha256` there, then:

```bash
rsync -a $B/report $B/SOURCE.sha256 vast:/root/grm-incoming-rev7/
# on the host
bash $G/publish.sh /root/grm-incoming-rev7/report /root/grm-incoming-rev7/SOURCE.sha256 rev7-restored          # dry run
bash $G/publish.sh /root/grm-incoming-rev7/report /root/grm-incoming-rev7/SOURCE.sha256 rev7-restored --apply
```

`publish.sh` verifies the source tree and checks its index against `GRM_EXPECT_*` first, refuses an existing name, copies into `published/NAME.partial`, verifies the copy, makes it read-only and renames it into place. For a report with another dataset count, set `GRM_EXPECT_DATASETS` (for example 275) for `publish.sh` and `stage.sh`.

## Build and stage a release

On the host that will serve it (`.cargo/config.toml` builds for the native CPU):

```bash
export GRM_EXTRA_PATH=/root/.cargo/bin:/root/grm/.tools/node/bin
bash $G/build-release.sh ui-r3-abc1234 abc1234           # dry run: commit, new paths, free space
bash $G/build-release.sh ui-r3-abc1234 abc1234 --apply
```

`GRM_SOURCE_REPO` selects the git checkout to archive from (default: the checkout holding the scripts, so set it when the scripts were copied out with `git archive`; `/root/grm` is a checkout on the host), `GRM_BUILD_JOBS` the cargo jobs (8) and `GRM_MIN_FREE_GIB` the free space required (40). The build writes `src-LABEL`, a temporary `target-LABEL` and `release/LABEL`, all of which must be new. Creating `src-LABEL` reserves the label before anything is extracted, so a second build under the same label refuses at once instead of sharing the source or target directory. A failed build leaves its paths for inspection, so retry under a new label.

Stage the candidate on a spare port. The browser check runs `web/scripts/browser-check.cjs` from the release's source with Playwright's Chromium:

```bash
export GRM_PLAYWRIGHT_CORE=/root/grm/.tools/browser-check/node_modules/playwright-core
export LD_LIBRARY_PATH=/root/grm/.tools/browser-libs/usr/lib/x86_64-linux-gnu:/root/grm/.tools/browser-libs/usr/lib/x86_64-linux-gnu/gbm
bash $G/stage.sh ui-r3-abc1234 rev7-20261025 8650           # dry run
bash $G/stage.sh ui-r3-abc1234 rev7-20261025 8650 --apply
```

It refuses the live port or a busy one, serves the pair on 127.0.0.1 as a [supervised start](#supervised-starts), requires the started process to own the port, runs `smoke.sh` against `GRM_EXPECT_*` and the browser check, stops the candidate, and keeps receipts under `stage/`. The candidate is stopped on every exit, including an error or INT, TERM or HUP at any point of the start; a candidate that ignores SIGTERM is killed and fails the stage, and if it cannot be killed the script says so and exits non-zero. If a stage is killed outright, its candidate keeps running and `serve-STAGE_PORT.starting` names it: the dry run reports it, and the next `stage.sh --apply` on that port stops it before staging. Only a full pass writes `go/LABEL--NAME`, atomically. `GRM_SKIP_BROWSER=1` approves on the HTTP checks alone and the approval records `browser=skipped`.

## Activate

```bash
bash $G/activate.sh ui-r3-abc1234 rev7-20261025             # dry run
bash $G/activate.sh ui-r3-abc1234 rev7-20261025 --apply
```

Activation refuses unless the candidate is approved and intact, exactly one process listens on the port, and that process is a release of this root serving one of its publications, from the root, and that pair passes the same checks, approval included, because it is the rollback target. It also refuses while `serve-PORT.starting` records an unresolved interrupted start, and when `serve-PORT.starting`, `serve-PORT.pid`, `logs/` or `logs/serve-PORT.log` is a link, directory or special file. A failed check exits 1 and a live service or path that is not as required exits 2, in both cases before anything changes. Under a lock it then:

1. records the live pair in `state/previous` and keeps it in `state/active`, each written and read back as described in [Approvals and state](#approvals-and-state), so that if the script is killed outright, `restore.sh --apply` brings the live pair back; if either write fails, nothing stops (exit 2);
2. stops the old process; if it does not stop within 10 s and still listens, it keeps serving and nothing starts (exit 3);
3. starts the candidate as a [supervised start](#supervised-starts) and requires the candidate's own PID to own the port and pass `smoke.sh`;
4. records the candidate in `state/active`, then its PID in `serve-PORT.pid`, then removes the start journal.

From step 2 until the script returns, any failure, script error (such as a failed write) or INT, TERM or HUP, at any point including between the fork and the assignment of the candidate's PID, stops every candidate process and restarts the old pair, which must in turn own the port under its new PID and pass `smoke.sh`; then the old pair serves again, `state/active` and `serve-PORT.pid` record it (exit 1). If that fails too, nothing healthy is serving (exit 4) and `deploy.log` has the details. A failed `serve-PORT.pid` write after step 4's state commit is a warning: the committed pair keeps serving (exit 0). The server never inherits the lock, so the next activation or rollback can take it.

### Supervised starts

`activate.sh`, `restore.sh`, `stage.sh` and the rollback restart in `activate.sh` start a server the same way (`start_service` in `common.sh`, built on `lib.sh`):

1. A new marker `spawn/TOKEN.XXXXXX` is created and the script takes a `flock` on it. `serve-PORT.starting` is written atomically with the fresh random token, the marker's name and device and inode, the starting script's PID and start time, and the pair. If either cannot be written, nothing starts.
2. The child is forked with the marker open. From that instant it, the `grm` it execs, and anything they fork hold the open file description that carries the lock, so the marker identifies every process of the start by device and inode (even if the marker is renamed or deleted), and the lock proves when none is left, including processes the script cannot see.
3. The script rewrites the journal with the child's PID and start time. The child waits for that registration and checks that its script is still alive, again just before exec; a child whose script died before registering, or that is resumed after its script died, exits without exec'ing. `grm` runs under `setsid`, so its process group is stopped with it, and with the token in its environment (`GRM_DEPLOY_SPAWN`).

An interruption at any point, even between the fork and the registration, stops every holder of the marker and the registered child (TERM, with CONT so a stopped process acts on it, then KILL after 10 s), and the marker is removed once its lock is free. If the script is killed outright, the journal stays and names a dead owner: `activate.sh` and `rollback.sh` refuse until it is resolved, and `restore.sh --apply` (live port) or `stage.sh --apply` (its port) resolves it. If `serve-PORT.pid` records the journal's token for a live PID that holds the marker, the start had committed and only the journal goes. Otherwise it stops every holder and the registered child, and removes the journal and marker only when the marker's lock is free; while a process it cannot see holds the lock, or the marker is gone, it keeps the journal, starts nothing and fails (exit 1). While the journal's owner is alive, `activate.sh`, `restore.sh` and `stage.sh` refuse to act on that port. `status.sh` reports the journal and the processes holding its marker.

The tests in [Checks](#checks) send TERM, INT, HUP, KILL or STOP through hooks (`GRM_TEST_SIGNAL_AT` and `GRM_TEST_SIGNAL` in the script, `GRM_TEST_CHILD_AT` in the forked child; inert unless set) at each step: after the old service stops, after the journal is written, after the fork, after the registration, after the PID is assigned, after the health check, after the state commit, after `serve-PORT.pid` is written, during a record write and after its rename, during the rollback restart, and in the child at the fork and just before exec.

After activating, check the public site from a browser or with the hosting browser check; Cloudflare's Browser Integrity Check rejects some non-browser clients.

## Roll back

```bash
bash $G/rollback.sh             # dry run
bash $G/rollback.sh --apply
```

This activates the pair in `state/previous` with every activation check (`activate.sh --rollback`), so a second rollback switches forward again. It needs a running managed service.

The target must verify; the live pair need not. If the live pair fails its checks (a damaged current release, a deleted or modified `dist/`, or a deleted or replaced running executable), the rollback still replaces it, but only when its recorded start identifies it, and it is exactly the pair `state/active` records: `serve-PORT.pid` must name that PID with its start time, the token in its environment and a marker under `spawn/` that still exists and that the process holds open, and the process must be `grm`, run from the root and serve the recorded publication. The executable's path is not used for this, so a server whose executable was deleted is still identified, and a process that merely runs from the root with the right arguments and a forged token or record is not. Otherwise it refuses with nothing stopped (exit 2). A server started before the spawn markers existed has no marker in `serve-PORT.pid`, so a damaged one is refused; replace it by hand. In that case the damaged pair is not written to `state/previous` and is never restarted: if the target then fails, nothing serves the port (exit 4), and `state/active` still names the damaged pair, which an argument-free `restore.sh` refuses. Bring the site back with an explicit `restore.sh LABEL PUBLICATION --apply` of an intact approved pair, such as the rollback target. The same applies when such a rollback is killed outright after the stop. While the rollback holds the lock, any explicit `restore.sh --apply` is refused.

To go back further, or before `state/previous` exists, activate an older approved pair directly; releases and publications are never deleted by these scripts:

```bash
bash $G/activate.sh ui-r2-5a2e522 rev7-20261025 --apply
```

The legacy `/root/grm-deploy-20261025/rollback.sh` on the host restores the revision-6 binary from `/root/grm` with its revision-6 data. Use it only to return to revision 6, not to undo a dashboard release.

## After a container restart

A restart stops `grm serve`. `/root/onstart.sh` starts only the tunnel connector, so the public site returns an origin error until the app is restored:

```bash
bash $G/restore.sh ui-r2-9f193bf rev7-20261025               # dry run; first time on a root without state/active
bash $G/restore.sh ui-r2-9f193bf rev7-20261025 --apply
bash $G/restore.sh --apply                                  # later: restores state/active
```

`restore.sh` verifies the pair and its approval (without arguments, also that `state/active` matches them), and refuses when another process holds the port. It first resolves an interrupted start on the port (see [Supervised starts](#supervised-starts)); the dry run only reports it. When the pair already serves and passes `smoke.sh`, it restarts nothing; if `state/active` or `serve-PORT.pid` does not record it yet, `--apply` records it, so a later argument-free restore finds it. Otherwise it refuses, before starting anything, when `serve-PORT.pid`, `serve-PORT.starting`, `state/active`, `logs/` or the server log is a link, directory or special file; then starts the pair as a supervised start, requires the started PID to own the port and pass `smoke.sh`, and records `state/active` and then `serve-PORT.pid`. A failed or interrupted restore stops what it started, leaving nothing on the port and `state/active` as it was (exit 1). A failed `serve-PORT.pid` write after the state commit is a warning: the healthy recorded pair keeps serving (exit 0), and the next `restore.sh --apply` records the PID. Because it is idempotent, a boot hook may run it before the connector; adding that to `/root/onstart.sh` is a manual decision:

```bash
GRM_DEPLOY_ROOT=/root/grm-deploy-20261025 bash /path/to/scripts/deployment/restore.sh --apply >>/var/log/grm-restore.log 2>&1
```

Then confirm that the connector is running (`pgrep -af cloudflared`) and that the public site loads.

## Tunnel connector

The tunnel token is a secret. Keep it only in a file owned by root with mode 0600, inside a 0700 directory, and never print, copy or commit it. `connector.sh` never reads it; cloudflared reads it through `--token-file`.

```bash
export GRM_CONNECTOR_DIR=/root/grm-cloudflare-20261026          # cloudflared, cloudflared.sha256, logs/
export GRM_TUNNEL_TOKEN_FILE=$GRM_CONNECTOR_DIR/tunnel.token     # the default
bash $G/connector.sh check                                       # binary hash and token-file permissions; changes nothing
nohup setsid bash $G/connector.sh run >>$GRM_CONNECTOR_DIR/supervisor.log 2>&1 </dev/null &   # starts the tunnel now
bash $G/connector.sh stop                                        # stops it now
```

`cloudflared.sha256` must hold exactly one complete line, `SHA256  cloudflared` (or `SHA256 *./cloudflared`) ending in LF, for a regular executable file (not a link) named `cloudflared`, and no other byte: another path, a second or blank line, bytes after the newline, a missing final newline, a CRLF ending, a NUL or other control byte fails. `supervisor.pid` must likewise be one complete line of a PID and a start time. `run` checks the binary again before every start, supervises cloudflared with `--protocol http2` (`GRM_TUNNEL_PROTOCOL`), restarting it 5 s after it exits, and writes `supervisor.pid`, `connector.pid`, `connector.lock`, `connector.spawn` and `logs/` in the connector directory. It refuses, without opening anything, when any of those, or an entry of `logs/`, is a link, FIFO or other special file. It holds `connector.lock` itself (cloudflared does not inherit it), so it cannot run beside the existing `/root/grm-cloudflare-20261026/run.sh` supervisor, which `/root/onstart.sh` starts today and `stop.py` stops.

cloudflared is supervised like a server start: the supervisor locks `connector.spawn` and forks with it open, so cloudflared holds that lock from the fork on, and runs under `setsid`. The forked child execs only while its supervisor is alive. TERM, INT or HUP at any point, including between the fork and the assignment of the child's PID (the handler uses bash's `$!`, which the fork itself sets), stops the child with its process group and exits 0; the supervisor never waits on a child it does not track. If the supervisor is killed outright, cloudflared keeps running and holding `connector.spawn`: the next `run` stops it before starting its own, and `stop` stops it after reporting the stale record, both only while no supervisor holds `connector.lock`; while the lock shows a holder they cannot see, `run` starts nothing (exit 1).

`stop` signals the PID in `supervisor.pid` only if it still has the recorded start time, works in this connector directory, and one of its descriptors is the open file description holding the flock on this directory's `connector.lock` (the kernel reports it in `/proc/PID/fdinfo`); the lock must be held, which `stop` checks with a non-blocking `flock`. A command line naming this script, an open descriptor on `connector.lock` without the lock, a lock on another directory's file, or a stale or reused PID gets no signal. Switching the boot hook from `run.sh` to `connector.sh run` is optional and manual. Stopping the connector takes the public site offline but leaves `grm serve` running.

To rotate the token, issue a new one for the tunnel in Cloudflare, write it under `umask 077` to a new file in the connector directory without echoing it or putting it on a command line, `mv` it over the token file, run `connector.sh check`, restart the connector, and revoke the old token in Cloudflare.

## Checks

```bash
printf '%s\n' scripts/deployment/*.sh scripts/deployment/tests/*.sh | xargs -n1 bash -n
bash scripts/deployment/tests/dry-run.sh
bash scripts/deployment/tests/apply.sh
```

`bash -n a.sh b.sh` parses only `a.sh` (the rest are its arguments), so the first line checks each file separately.

Both suites delete their temporary root on exit, so they take it from `make_test_root` (`tests/tmproot.sh`), which stops the suite (exit 2) when `mktemp` fails and refuses a directory that is `/`-level, not empty, or contains the repository. The earlier `T=$(cd -- "$(mktemp -d)" && pwd -P)` fell back to the current directory when `mktemp` failed, because `cd -- ""` succeeds; run from the repository with an unusable `TMPDIR`, the exit trap deleted the whole working tree. From Git Bash on Windows, set `MSYS_NO_PATHCONV=1` before passing Linux paths through `wsl.exe`: without it Git Bash rewrites `TMPDIR=/home/...` into `TMPDIR=C:/Program Files/Git/home/...`, which is how that `TMPDIR` became unusable.

Both suites count their assertions and fail when the count differs from the expected total (`EXPECTED_ASSERTIONS`), and an EXIT trap prints a `FAIL` line when a command outside an assertion ends the suite early, so a SIGPIPE or other failure under `set -e` and `pipefail` can no longer skip the rest silently. Pipelines in the suites read their whole input (`sed -n`, not `head` or `grep -q` on a pipe); `tests/dry-run.sh` checks that on a 2.7 MB pipeline, where `head` is killed by SIGPIPE. Every guard against a hang is `timeout -k 5 N`: a script that traps TERM and is blocked opening a FIFO never runs its trap, so a TERM alone would leave it, and the suite, waiting forever. The mutation run found that gap: with the connector FIFO defect reintroduced (case 31), the suite hung instead of failing.

`tests/dry-run.sh` builds a throwaway root in a temporary directory with a two-dataset report, stub releases and a stub cloudflared. It requires ports 1 and 2 to have no listener and opens none, never runs the stub `grm`, `cargo` or `npm`, and runs only `publish.sh --apply` into the temporary root. It checks the refusals: damaged, unlisted or missing files; symlinks, FIFOs, linked directories and stray manifests in releases and reports; wrong schema, revision or count; missing, legacy, tampered, repeated, unknown or mismatched approval fields, including another report copy and names that differ only where a regex would match; approvals and legacy `GO-LABEL` files with a NUL followed by a conflicting stage, commit, publication, report or release manifest, schema, revision or dataset count, and with CR, tab, DEL or non-ASCII bytes, a missing final newline, blank lines or stray spaces, or that are empty, FIFOs or directories; state files with a NUL, blank line, CRLF, repeated field, missing final newline or missing fields, or that are FIFOs, directories or links; start journals that are malformed, directories or owned by a live process, and FIFOs at `serve-PORT.pid` or the server log; the live port; path-like names; a readable token file; a changed binary; a cloudflared manifest naming another, path-like or absolute path, with an extra complete line, unterminated junk, no final newline, a NUL, CRLF, a blank line, a tab, an uppercase hash or nothing in it (the valid line and `*./cloudflared` pass); connector `stop` with malformed, two-line, NUL, FIFO, stale, reused or foreign PIDs, an unlocked or FIFO `connector.lock`, and a FIFO manifest; and start journals without a spawn marker, with a marker of another token, a half child registration or an unknown field. It also checks that dry runs, with and without a start journal, leave the root unchanged, that an interrupted start is reported and blocks activation and rollback, and that the connector check never prints the token.

`tests/apply.sh` needs `cc` and two free loopback ports, which it picks, and takes about 15 to 20 minutes on a loaded 12-core host. Started where SIGINT or SIGHUP is ignored (a background job, `nohup`), it restarts itself with the default dispositions, because a shell cannot trap a signal that was ignored when it started. It compiles a stub `grm` that serves the app and report on 127.0.0.1 and, by environment, exits at once, listens late, serves a changed page or ignores SIGTERM. It runs `publish`, `stage`, `restore`, `activate`, `rollback`, `approve-legacy`, `status`, `build-release` and `connector run`/`stop` with `--apply` in a temporary root: 237 to 275 to 237 to 275 datasets and back after a simulated crash, with a contradictory `GRM_EXPECT_DATASETS` set; that the server never holds `deploy.lock`; refusal when the live pair lacks a modern approval, and the legacy migration; before the stop, refusal with the old service untouched when `state/previous` or `state/active` is a directory (nothing is moved into it), a link, a dangling link (its target is not created), a FIFO, in an unwritable `state/` (skipped as root), truncated in its temporary file, or changed after its rename, and when `serve-PORT.starting`, `serve-PORT.pid` (directory, link, FIFO), the server log or `deploy.lock` is a special file; a FIFO `deploy.log` that is skipped; a write killed before its rename that leaves the old record; TERM at each of seven steps of activation (from after the old service stops to after `serve-PORT.pid` is written), INT before and HUP after the PID assignment, each restoring the old pair with no candidate process or journal left; KILL at each step, after which an argument-free restore stops the orphaned candidate by its token and restores `state/active`, or keeps a start that had committed; KILL during the rollback restart; a state commit and a `serve-PORT.pid` write that fail after the stop; TERM at six steps of a restore, leaving nothing on the port, and KILL at four, resolved by the next restore; `serve-PORT.pid` as a directory, link or FIFO refused before a restore starts anything, a read-only one replaced, and one that fails at its exact boundary, after the state commit, keeping the healthy pair; a failed health check and an external SIGTERM after the stop, each restoring the old pair; an external SIGKILL after the stop followed by an argument-free restore; an impostor process serving identical bytes, which fails both the candidate and the rollback because neither started PID owns the port; a rollback from a damaged current release to the intact previous pair, refused when `state/active` names another pair, refusing an explicit restore while it holds the lock, never returning to the damaged pair, leaving nothing serving when the target fails and when killed, then refusing an impostor, refused by an argument-free restore and recovered by an explicit one; restore recording a healthy service that `state/active` lacks; staging with a candidate that ignores SIGTERM, with an interrupted stage, with TERM at four steps and KILL at two (the dry run reports the orphan, `--apply` stops it), with receipt directories for the current second already present, and with an approval path that is or becomes a directory, leaving no process or listener; two concurrent builds of one label; and connector supervisors that keep the lock from cloudflared, ignore a foreign or reused PID, and stop on a changed binary. Since the third review it also covers: TERM at the registration step of activation, restore and stage; KILL right after the fork (the unregistered child exits without exec'ing, and the journal is resolved) and after the registration; a child paused at the fork or just before exec while activation, rollback, restore or stage is killed outright, then either resolved (the resolving script stops it, and resuming it afterwards starts nothing) or resumed first (it exits without exec'ing and the journal is kept until resolved); a non-dumpable process holding an interrupted start's marker lock, and a deleted marker, both of which keep the journal and start nothing; a server whose running executable was deleted, replaced by rollback by its recorded identity, with a record naming another start time, PID, token, marker or pair, no record, a missing marker, an impostor that does not hold the marker, a plain activation and a restore all refused with nothing stopped, and a server whose `dist/` was deleted; and for the connector, TERM at each of five supervisor steps and INT and HUP between the fork and the PID assignment (each exits 0 with nothing left), a child paused before exec while its supervisor is killed (resumed, or stopped by `stop`), a killed supervisor's cloudflared stopped by the next `stop` and by the next `run`, a hidden holder of `connector.spawn` that stops `run` from starting anything, FIFOs at `connector.lock`, `connector.spawn`, `supervisor.pid`, `connector.pid` and in `logs/` and links at `supervisor.pid` and `logs/` refused without blocking, and forged supervisors (a shell in the directory with this script's path and `connector.lock` open on descriptor 9, with the lock free or held by the real supervisor; a holder of another directory's lock; the lock holder working in another directory; the real PID with another start time) that get no signal. It kills every process it started on exit, including non-dumpable ones, which it tracks by PID.

`data/work/deployment-mutation-20261027/grm-mutate-v3.py` reintroduces each of 33 defects (the 32 found by three reviews, plus a stage receipt-directory collision found during the work) into a temporary copy of `scripts/deployment/` and runs the matching suite; a defect counts as caught only when the suite exits non-zero with an explicit `FAIL` line and does not time out. Its README lists the runs and their receipts; the earlier harnesses and their receipts are kept unchanged.
