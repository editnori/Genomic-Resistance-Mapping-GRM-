// A tab opened before a switch to a new app label still asks for the old app's
// hashed chunks under /assets/. The Worker serves those from state/previous's
// app, and only those: other paths, /data/ and a same-label switch are as
// before, and a state/previous that doesn't parse leaves the miss to
// index.html. No account and no Miniflare.
// usage: node ops/cloudflare/tests/previous-assets.mjs
import assert from 'node:assert/strict';
import { handle } from '../worker/src/index.js';

const INDEX = '<!doctype html><title>b</title>';
const OLD = 'app/ui-a/dist/assets/PresentMode-old.js';
const base = {
  'app/ui-a/dist/index.html': '<!doctype html><title>a</title>',
  [OLD]: 'export const old = 1;',
  'app/ui-a/dist/favicon-a.svg': '<svg/>',
  'app/ui-b/dist/index.html': INDEX,
  'app/ui-b/dist/assets/PresentMode-new.js': 'export const now = 1;',
  'pub/pub-a/report/only-a.json': '{"a":1}',
  'pub/pub-b/report/index.json': '{"schema":2}',
};
const ACTIVE = JSON.stringify({ app: 'ui-b', publication: 'pub-b' });
// state/previous as stored text; undefined leaves it out
const bucketOf = previous => {
  const objects = { ...base, 'state/active': ACTIVE };
  if (previous !== undefined) objects['state/previous'] = previous;
  const lookup = key => {
    const value = objects[key];
    if (value === undefined) return null;
    return { size: Buffer.byteLength(value), body: value, httpEtag: `"${key}"`, text: async () => value };
  };
  return { get: async key => lookup(key), head: async key => lookup(key) };
};
const pair = (app, publication) => JSON.stringify({ app, publication });

const switched = bucketOf(pair('ui-a', 'pub-a'));
const appOnly = bucketOf(pair('ui-a', 'pub-b'));
const sameLabel = bucketOf(pair('ui-b', 'pub-a'));
const first = bucketOf(undefined);
const missingApp = bucketOf(pair('ui-gone', 'pub-a'));
const corrupt = ['not json', 'null', '{}', pair('../x', 'pub-a'), JSON.stringify({ app: 'ui-a' })].map(bucketOf);

const IMMUTABLE = 'public, max-age=31536000, immutable';
const asIndex = { status: 200, type: 'text/html', cache: 'no-cache', body: INDEX };
// [bucket, method, path, request headers, expected]
const cases = [
  [switched, 'GET', '/assets/PresentMode-new.js', {}, { status: 200, type: 'text/javascript', cache: IMMUTABLE, body: 'export const now = 1;' }],
  [switched, 'GET', '/assets/PresentMode-old.js', {}, { status: 200, type: 'text/javascript', cache: IMMUTABLE, body: 'export const old = 1;' }],
  [switched, 'HEAD', '/assets/PresentMode-old.js', {}, { status: 200, type: 'text/javascript', cache: IMMUTABLE }],
  [switched, 'GET', '/assets/PresentMode-old.js', { Range: 'bytes=0-5' }, { status: 206, type: 'text/javascript', range: 'bytes 0-5/21' }],
  [switched, 'GET', '/assets/PresentMode-old.js', { 'If-None-Match': `"${OLD}"` }, { status: 304, cache: IMMUTABLE }],
  [appOnly, 'GET', '/assets/PresentMode-old.js', {}, { status: 200, type: 'text/javascript', cache: IMMUTABLE, body: 'export const old = 1;' }],
  [switched, 'GET', '/assets/missing-everywhere.js', {}, asIndex],
  [switched, 'GET', '/favicon-a.svg', {}, asIndex],
  [switched, 'GET', '/data/only-a.json', {}, { status: 404, type: 'text/plain', cache: 'no-store' }],
  [sameLabel, 'GET', '/assets/PresentMode-old.js', {}, asIndex],
  [first, 'GET', '/assets/PresentMode-old.js', {}, asIndex],
  [missingApp, 'GET', '/assets/PresentMode-old.js', {}, asIndex],
  ...corrupt.map(bucket => [bucket, 'GET', '/assets/PresentMode-old.js', {}, asIndex]),
  ...corrupt.map(bucket => [bucket, 'GET', '/assets/PresentMode-new.js', {}, { status: 200, type: 'text/javascript', cache: IMMUTABLE }]),
];
const quiet = console.error;
console.error = () => {};  // the corrupt cases log their parse errors
try {
  for (const [bucket, method, path, headers, want] of cases) {
    const response = await handle(new Request('https://grm.test' + path, { method, headers }), bucket);
    const what = `${method} ${path} ${JSON.stringify(headers)}`;
    assert.equal(response.status, want.status, what);
    if (want.type) assert.ok(response.headers.get('Content-Type').startsWith(want.type), `${what}: ${response.headers.get('Content-Type')}`);
    if (want.cache) assert.equal(response.headers.get('Cache-Control'), want.cache, what);
    if (want.range) assert.equal(response.headers.get('Content-Range'), want.range, what);
    if (want.body !== undefined) assert.equal(await response.text(), want.body, what);
  }
} finally {
  console.error = quiet;
}
console.log(`previous-assets: ${cases.length}/${cases.length}`);
