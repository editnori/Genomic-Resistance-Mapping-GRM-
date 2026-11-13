# Static hosting on Cloudflare R2

From revision 9 on, genomicresistance.com is meant to be served with no origin process: a Cloudflare Worker (`ops/cloudflare/worker/`) answers every request from one R2 bucket, `grm-site`. The bucket holds immutable app releases and report publications, and one pointer object, `state/active`, names the pair the Worker serves. The release scripts in `ops/cloudflare/scripts/` write the bucket. This replaces the `grm serve` host described in [deployment.md](deployment.md), whose toolkit stays for a host if one is ever needed again.

The Worker serves genomicresistance.com this way now. On 13 November 2026 it answered with the revision-14 report, 274 datasets over 23 species. The tests under [Checks](#checks) ran locally with Miniflare (`wrangler dev --local`) on synthetic trees.

## What the Worker serves

`ops/cloudflare/worker/src/index.js` matches `grm serve` (`src/serve.rs`) route by route:

| Request | Answer |
|---|---|
| `/data/PATH` | `pub/NAME/report/PATH` of the active publication; 404 `not found` (`no-store`) when it is missing, a directory, or the path has a `..` segment. `/data/research.json` is 404, as the app expects |
| `/` and any other path | `app/LABEL/dist/PATH` of the active app; an unknown path gets `index.html`, so client routes such as `/species/NAME` load. A path under `/assets/` that the active app lacks is tried in the `state/previous` app first, when that app has another label |
| Methods other than GET and HEAD | 405 with `Allow: GET, HEAD` |
| `Range: bytes=…` | one range: 206 with `Content-Range`; a list, a malformed or unsatisfiable range: 416 with `Content-Range: bytes */SIZE`. Ranges past the end are clamped, as in `grm serve` |
| Content types | `serve.rs`'s table (html, js, mjs, css, json, svg, png, woff2, woff, ico), else `application/octet-stream` |
| `Cache-Control` | `public, max-age=31536000, immutable` for files under `/assets/` that exist, in the active app or the previous one; `no-cache` for everything else served, including the `index.html` fallback; `no-store` for errors |

Beyond `grm serve`, every object response carries an `ETag` and honours `If-None-Match` (304) and `If-Range`, and the Worker answers 503 while no pair is active. One difference cannot be removed: the Workers runtime resolves `..` and `%2e%2e` segments in the URL before the Worker sees it, so `/data/../x` arrives as `/x` and gets `index.html` (200) where `grm serve` answers 404. Nothing outside the active pair's prefixes, and the previous app's `assets/`, is reachable either way. Browsers resolve those segments themselves, so the app never sends them.

Each request reads `state/active` once, so it is answered from that one pair; only an `/assets/` miss also reads `state/previous` (below). R2 is strongly consistent, so a switch takes effect for the next request everywhere.

**A tab opened before a switch.** It keeps the old app in memory and asks for that app's lazy-loaded chunks under `/assets/`, such as the Present page's. The new app doesn't have them, so the Worker serves them from the `state/previous` app (hashed names can't collide). That covers one activation back, since `state/previous` names the pair that was active just before the latest activation. After an app switch and then a publication-only switch, both records name the new app, so a tab from before the app switch gets `index.html` for those chunks, as does a tab two app switches old, and its Present page fails until it reloads. A `state/previous` that doesn't parse is ignored here, so those misses get `index.html` rather than an error. Until it reloads, such a tab also reads `/data/` from the new publication under the index it already holds.

## Bucket layout and guarantees

| Key | Contents |
|---|---|
| `app/LABEL/CLAIM` | reserves the label; written first |
| `app/LABEL/dist/…`, `app/LABEL/COMMIT` | the built `web/dist`, and the commit it came from |
| `app/LABEL/MANIFEST.sha256` | `SHA256  ./path` for `dist/`; written last, the seal |
| `pub/NAME/CLAIM` | reserves the name, with the schema, revision and dataset count checked; written first |
| `pub/NAME/report/…` | the report tree |
| `pub/NAME/SOURCE.sha256` | its manifest; written last, the seal |
| `state/active`, `state/previous` | one JSON line each: the served pair, and the rollback target |
| `state/history/TIME--LABEL--NAME` | one record per activation, never rewritten |

- **A release is immutable once published.** `publish.sh` and `publish-app.sh` refuse a name whose claim or seal exists, so an interrupted upload leaves its claim and the retry takes a new name. Nothing in these scripts deletes or overwrites under `app/` or `pub/`. Two operators publishing the same new name at the same instant are not guarded against (Wrangler has no conditional write); run one operator at a time. The same holds for `activate.sh` and `rollback.sh`: two at once could each record the other's pair as `state/previous`.
- **Only sealed, intact releases are activated.** `activate.sh` requires both seals, requires each claim to name its seal's digest, and reads back and hashes every object of both before it switches.
- **Activation is atomic.** It writes `state/previous` (the old active pair), a history record, then `state/active`, each read back and compared. The Worker reads `state/active`, one object, so each response comes from one pair, the old or the new (an `/assets/` miss also reads `state/previous`, above). With `GRM_SITE_URL` set, `smoke.py` then checks the served site, and a failure restores both state objects byte for byte (exit 1).
- **Rollback names the previous release.** `rollback.sh` activates the pair in `state/previous`, with every activation check, and swaps the two, so a second rollback switches forward again.
- **Every step is logged.** Each `--apply` step appends a UTC line to `GRM_SITE_LOG`, and failures too. Dry runs (no `--apply`) run every read-only check, print what they would do and log nothing.

## Scripts

All take their configuration from the environment; see the header of `ops/cloudflare/scripts/site-lib.sh`.

| Script | Does |
|---|---|
| `publish.sh REPORT MANIFEST NAME [--apply]` | runs `scripts/deployment/publish.sh`'s checks as a dry run on an empty root (tree against manifest, index against `GRM_EXPECT_*`), then claim, upload, read back every object, seal. `GRM_EXPECT_REVISION` and `GRM_EXPECT_DATASETS` are both required |
| `publish-app.sh DIST LABEL COMMIT [--apply]` | the same for a built `web/dist` |
| `activate.sh LABEL NAME [--apply]` | the switch described above |
| `rollback.sh [--apply]` | activate `state/previous` |
| `status.sh` | prints `state/active`, `state/previous` and the last log lines |
| `smoke.py URL APP_MANIFEST REPORT_MANIFEST SCHEMA REVISION DATASETS [--sweep]` | the served check; `--sweep` fetches every report file |
| `wrangler.sh ARGS…` | the pinned Wrangler with the token, for bucket creation and deploy |

## Checks

`ops/cloudflare/tests/`, run against Miniflare state with no account:

- `release.sh`: publish refusals (missing expectations, wrong count or revision, tampered source, bad name, existing name, existing claim), activation refusals (unsealed, tampered object, missing app, already active), activation and rollback state, the served check and its restore, the log.
- `routes.py`: the Worker against `grm serve` on the same tree, route by route: status, the headers above and body SHA-256.
- `ranges.py`: range requests on one large file against its bytes on disk.
- `long-key.mjs` and `previous-assets.mjs` (plain `node`, a stub bucket): a path whose key passes R2's 1,024-byte limit answers as missing, not 500; and an `/assets/` miss is served from the previous app only when that app has another label (with ranges and conditional requests as usual), never for `/data/` or other paths, and a `state/previous` that doesn't parse gives `index.html`, not 500.

## Runbook: revision 9 live

Roles as in the recovery runbook: Layth holds the Cloudflare account and the token; main runs the commands; nothing is printed that holds the token.

### 1. The token

Create one API token in the dashboard (My Profile › API Tokens › Create Custom Token), **restricted to the one account and the one zone**, with a short expiry (for example 7 days) and, if possible, Client IP Address Filtering to the operator's address:

| Scope | Permission | Why |
|---|---|---|
| Account › Workers R2 Storage | Edit | create `grm-site`, write releases and `state/` objects. Cloudflare API tokens cannot narrow this to one bucket; bucket-scoped R2 tokens exist only for the S3 API, which these scripts do not use |
| Account › Workers Scripts | Edit | deploy the `grm-site` Worker |
| Zone `genomicresistance.com` › Workers Routes | Edit | attach the Worker's custom domain |
| Zone `genomicresistance.com` › DNS | Edit | remove the tunnel's CNAME, and let the custom domain create its record |

If `wrangler deploy` refuses the custom domain for want of a zone read, add Zone › Zone › Read on the same zone; nothing else should be needed, since `CLOUDFLARE_ACCOUNT_ID` is always set. For later releases only, a second token with Workers R2 Storage › Edit alone is enough; revoke the first once the site is live.

**Passing it.** Main asks Layth for it through a private secret card with file delivery. The file must be mode 0600 and owned by the operator; `site-lib.sh` checks both and refuses anything else. If the delivered file is not 0600, copy it once with `umask 077; cat "$DELIVERED" >"$HOME/.grm-site-token"` and point `GRM_SITE_TOKEN_FILE` there. The scripts read the file into Wrangler's environment only, never onto a command line, into a log or to the terminal. `site-lib.sh` turns shell tracing off, so even `bash -x` doesn't print it. When done: `shred -u "$GRM_SITE_TOKEN_FILE"` and revoke the token in the dashboard.

### 2. Setup

From a checkout of the commit that carries `ops/cloudflare/`, with Node 22 and Python 3:

```bash
npm --prefix ops/cloudflare/worker ci
export GRM_SITE_TARGET=remote GRM_SITE_ACCOUNT_ID=<account id> GRM_SITE_TOKEN_FILE=$HOME/.grm-site-token
export GRM_SITE_LOG=$PWD/site-$(date -u +%Y%m%d).log
W=ops/cloudflare/scripts
bash $W/wrangler.sh whoami                     # must name the GRM account
bash $W/wrangler.sh r2 bucket create grm-site
```

### 3. The app

Revision 9 goes live with main's app at `fd558ce2` (the decision of 5 Nov): its `web/` build is byte-identical to the one lane U-16 checked on revision 9 (1,669 checks), and REVIEW-HOST passed the current browser check through the Worker on revision 9 with it (1,694 checks). Build it from a clean archive:

```bash
git archive fd558ce262637c007d3a8224c3c6ef5bdb5cf11b web | tar -x -C /tmp/web-fd558ce2 && (cd /tmp/web-fd558ce2/web && npm ci && npm run build)
bash $W/publish-app.sh /tmp/web-fd558ce2/web/dist ui-main-fd558ce2 fd558ce262637c007d3a8224c3c6ef5bdb5cf11b            # dry run
bash $W/publish-app.sh /tmp/web-fd558ce2/web/dist ui-main-fd558ce2 fd558ce262637c007d3a8224c3c6ef5bdb5cf11b --apply
```

### 4. The publications

Revision 9 from `grm-release-rev9-275-20261030.tar.gz` (sha256 `a694744eda82678bc90afec1915b2da47522d73840af6f6c2c0948399e315ce5`), and revision 8 from `grm-release-rev8-273-20261029.tar.gz` as the rollback target. Each archive holds `SOURCE.sha256`, `SKIPS.json`, `RELEASE-NOTES.md` and `data/report/`.

```bash
sha256sum -c grm-release-rev9-275-20261030.tar.gz.sha256
mkdir /tmp/rev9 && tar -xzf grm-release-rev9-275-20261030.tar.gz -C /tmp/rev9
sha256sum /tmp/rev9/SOURCE.sha256              # must be 213a9978327bbbdf04bef2722937dc1b92b49729512778a35ee233544c7a677e
GRM_EXPECT_REVISION=9 GRM_EXPECT_DATASETS=275 bash $W/publish.sh /tmp/rev9/data/report /tmp/rev9/SOURCE.sha256 rev9-275-20261030           # dry run
GRM_EXPECT_REVISION=9 GRM_EXPECT_DATASETS=275 bash $W/publish.sh /tmp/rev9/data/report /tmp/rev9/SOURCE.sha256 rev9-275-20261030 --apply
# the same for revision 8 into /tmp/rev8: SOURCE.sha256 3aac8c60…a061, GRM_EXPECT_REVISION=8 GRM_EXPECT_DATASETS=273, name rev8-273-20261029
```

Each upload of about 1.1 GB, 300 objects, then reads every object back.

### 5. Activate, then deploy

Activate revision 8 first and then revision 9, so `state/previous` holds a real rollback target. The Worker is not deployed yet, so there is no served check here (leave `GRM_SITE_URL` unset):

```bash
bash $W/activate.sh ui-main-fd558ce2 rev8-273-20261029 --apply
bash $W/activate.sh ui-main-fd558ce2 rev9-275-20261030 --apply
bash $W/status.sh                              # active rev9, previous rev8
```

The custom domain cannot be created while `genomicresistance.com` has a CNAME, and the tunnel's route is one. In the dashboard (DNS › Records), Layth or main deletes the CNAME for `genomicresistance.com` that points at `<tunnel-uuid>.cfargotunnel.com` (and `www`, if one exists and should serve too). The site has answered 530 since 31 October, so this removes nothing that works. Then:

```bash
bash $W/wrangler.sh deploy                     # Worker grm-site, custom domain genomicresistance.com
```

Wrangler creates the DNS record and the certificate. `wrangler.toml` turns off `workers.dev` and preview URLs, so the Worker answers only on the domain.

Every hostname on this account must have a Cloudflare Access application, or Cloudflare answers 403 "error code: 1050" for it. Keep the existing application for `genomicresistance.com` (Bypass, Include Everyone).

`wrangler.toml` routes only the apex, `genomicresistance.com`, so no other hostname serves the site today. A new hostname (`www`, if it should serve) needs both of these before it is announced:
- **A route to the Worker.** Add it to `routes` in `ops/cloudflare/worker/wrangler.toml`, for example `{ pattern = "www.genomicresistance.com", custom_domain = true }`, and run `wrangler.sh deploy` again. Wrangler creates that hostname's DNS record and certificate, as it did for the apex.
- **Its own Access application.** In Zero Trust › Access › Applications, create a self-hosted application for that hostname with one policy: action Bypass, Include Everyone. The apex's application doesn't cover it.

### 6. Public check

```bash
(cd /tmp/web-fd558ce2/web/dist && find . -type f | LC_ALL=C sort | xargs sha256sum) >/tmp/app.sha256
python3 -I $W/smoke.py https://genomicresistance.com /tmp/app.sha256 /tmp/rev9/SOURCE.sha256 2 9 275 --sweep
```

Expect `smoke pass … schema=2 revision=9 datasets=275 species=23 skips_sha256=9f83f721…e746`, with every one of the 302 report files matching. Then the release's browser check against `https://genomicresistance.com` with the revision-9 tree (expect it to pass as at the revision-9 switch). Browser Integrity Check can challenge non-browser clients; if `smoke.py` gets 403 while a browser works, add a WAF skip rule for the operator's IP for the check, not a global change.

### 7. Rollback

```bash
GRM_SITE_URL=https://genomicresistance.com bash $W/rollback.sh             # dry run: must name rev8-273-20261029
GRM_SITE_URL=https://genomicresistance.com bash $W/rollback.sh --apply     # only on a failed public check
```

With `GRM_SITE_URL` set, a rollback that fails its own served check restores revision 9 and exits 1. Later activations (revision 11) run `activate.sh` with `GRM_SITE_URL` set from the start.

### 8. Retire the tunnel

- Zero Trust › Networks › Tunnels: delete the Vast-era tunnel. Its connector died with Vast, and deleting the tunnel invalidates its token.
- Zero Trust › Access › Applications: **keep** the application for `genomicresistance.com` (action Bypass, Include Everyone). It is not the tunnel's: every hostname on this account must have an Access application, and without one Cloudflare answers 403 "error code: 1050" for the whole site. Delete nothing here.
- DNS: confirm no record still points at `*.cfargotunnel.com`.
- Revoke the setup token and `shred -u` its file.

## Cost

Cloudflare's published free allowances (R2 pricing and Workers pricing pages, read 5 Nov 2026), for revisions 8 and 9 (about 1.1 GB each) plus an app of a few MB, about 2.3 GB:

| Item | Free allowance | This site |
|---|---|---|
| R2 Standard storage | 10 GB-month per month | about 2.3 GB-month; revision 11 adds about 1.1 GB |
| R2 Class A (PutObject, ListObjects) | 1 million per month | about 300 per publication, 5 per activation |
| R2 Class B (GetObject, HeadObject) | 10 million per month | 2 per plain GET (the pointer and the object), 3 per HEAD or range; about 330 per publication and per activation for the read-back |
| R2 egress | free | free |
| Workers Free requests | 100,000 per day, 10 ms CPU per invocation | 1 per request; one app visit makes a few dozen |

At the Workers Free ceiling every day, Class B use is about 6 to 9 million a month, inside the free 10 million. So the site costs nothing until it passes about 100,000 requests a day. Past that, Free returns error 1027 for the rest of the UTC day; Workers Paid ($5 a month, 10 million requests included, $0.30 per further million) removes the daily cap, and R2 past its allowances is $0.015 per GB-month, $4.50 per million Class A and $0.36 per million Class B. The Worker streams R2 bodies without touching them, so a 131.9 MiB report costs one request and little CPU; Cloudflare sets no response size limit.
