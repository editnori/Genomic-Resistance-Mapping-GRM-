// R2 keys are at most 1,024 bytes, so a path whose key would be longer names
// nothing: it answers as a missing file (404 under /data/, index.html for an
// app route), never 500. The bucket here throws on a long key, as R2 does.
// No account and no Miniflare. usage: node ops/cloudflare/tests/long-key.mjs
import assert from 'node:assert/strict';
import { handle } from '../worker/src/index.js';

const objects = new Map([
  ['state/active', JSON.stringify({ app: 'ui-t', publication: 'pub-t' })],
  ['app/ui-t/dist/index.html', '<!doctype html><title>t</title>'],
  ['pub/pub-t/report/index.json', '{"schema":2}'],
]);
const lookup = key => {
  if (new TextEncoder().encode(key).length > 1024) throw new Error(`R2 refuses a ${key.length}-byte key`);
  const value = objects.get(key);
  if (value === undefined) return null;
  return { size: Buffer.byteLength(value), body: value, httpEtag: '"t"', text: async () => value };
};
const bucket = { get: async key => lookup(key), head: async key => lookup(key) };

const long = 'a'.repeat(1100);
const cases = [
  ['GET', `/data/${long}.json`, 404],
  ['HEAD', `/data/${long}.json`, 404],
  ['GET', `/species/${long}`, 200, 'text/html'],
  ['GET', '/data/index.json', 200, 'application/json'],
  ['GET', '/species/x', 200, 'text/html'],
];
for (const [method, path, status, type] of cases) {
  const response = await handle(new Request('https://grm.test' + path, { method }), bucket);
  assert.equal(response.status, status, `${method} ${path.slice(0, 40)}…`);
  if (type) assert.ok(response.headers.get('Content-Type').startsWith(type), `${method} ${path.slice(0, 40)}…`);
}
console.log(`long-key: ${cases.length}/${cases.length}`);
