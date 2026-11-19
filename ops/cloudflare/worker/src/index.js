// Serves the GRM dashboard and the active report publication from one R2
// bucket, the way `grm serve` (src/serve.rs) serves them from disk:
//
//   /data/<path>  pub/<publication>/report/<path>; a missing path is 404
//   /<path>       app/<label>/dist/<path>; an unknown path gets index.html
//                 so client-side routes load. An /assets/ miss is tried in
//                 state/previous's app first, so a tab opened before an app
//                 switch still loads its lazy chunks
//
// state/active names the pair. ops/cloudflare/scripts/activate.sh writes it
// with one PUT, and every request reads it, so a switch is atomic: each
// response comes from the pair its request read, or, for an /assets/ miss,
// from state/previous's app. Releases under app/ and pub/ are never
// overwritten (publish.sh refuses an existing prefix).
//
// Beyond grm serve: an ETag on every object response, with If-None-Match
// (304) and If-Range, and 503 while no pair is active.

const ACTIVE_KEY = 'state/active';
const PREVIOUS_KEY = 'state/previous';
const NAME = /^[A-Za-z0-9][A-Za-z0-9._-]*$/;
const IMMUTABLE = 'public, max-age=31536000, immutable';

// serve.rs content_type, keyed by the last path component's extension.
const TYPES = {
  html: 'text/html; charset=utf-8',
  js: 'text/javascript; charset=utf-8',
  mjs: 'text/javascript; charset=utf-8',
  css: 'text/css; charset=utf-8',
  json: 'application/json',
  svg: 'image/svg+xml',
  png: 'image/png',
  woff2: 'font/woff2',
  woff: 'font/woff',
  ico: 'image/x-icon',
};

export function contentType(key) {
  const name = key.slice(key.lastIndexOf('/') + 1);
  const dot = name.lastIndexOf('.');
  return (dot > 0 && TYPES[name.slice(dot + 1)]) || 'application/octet-stream';
}

// serve.rs resolve: the path's components joined under `prefix`, or null when
// one would leave it. Like Rust's Path::components, empty segments and
// non-leading `.` drop out; a leading `.` or any `..` refuses the path. A
// trailing `/` names a directory, which is never a file there either.
export function keyUnder(prefix, urlPath) {
  if (urlPath.endsWith('/')) return null;
  const parts = urlPath.replace(/^\/+/, '').split('/');
  if (parts[0] === '.' || parts.includes('..')) return null;
  const kept = parts.filter(part => part !== '' && part !== '.');
  return kept.length ? prefix + kept.join('/') : null;
}

// serve.rs requested_range: one `bytes=` range, no lists. Repeated Range
// headers arrive joined with ", ", so they fail as a list does there.
export function requestedRange(header, size) {
  if (header === null) return { full: true };
  if (!header.startsWith('bytes=')) return null;
  const spec = header.slice(6);
  if (spec.includes(',')) return null;
  const dash = spec.indexOf('-');
  if (dash < 0) return null;
  const first = spec.slice(0, dash), last = spec.slice(dash + 1);
  if (first === '') {
    const suffix = digits(last);
    if (suffix === null || suffix === 0 || size === 0) return null;
    const length = Math.min(suffix, size);
    return { offset: size - length, length };
  }
  const start = digits(first);
  if (start === null || start >= size) return null;
  let end = size - 1;
  if (last !== '') {
    const parsed = digits(last);
    if (parsed === null) return null;
    end = Math.min(parsed, size - 1);
  }
  if (end < start) return null;
  return { offset: start, length: end - start + 1 };
}

// Rust's u64 parse: ASCII digits only, with an optional leading `+`, at most
// 2^64 - 1. Values past 2^53 lose precision, which no object size reaches.
function digits(text) {
  const match = /^\+?([0-9]+)$/.exec(text);
  if (!match || BigInt(match[1]) > 0xffffffffffffffffn) return null;
  return Number(match[1]);
}

function text(status, body, extra = {}) {
  return new Response(body, {
    status,
    headers: { 'Content-Type': 'text/plain; charset=UTF-8', 'Cache-Control': 'no-store', ...extra },
  });
}

// The request path as sent, before `?`: grm serve neither decodes nor
// normalises it, and keyUnder applies its rules.
function rawPath(url) {
  const afterScheme = url.indexOf('//') + 2;
  const slash = url.indexOf('/', afterScheme);
  const path = slash < 0 ? '/' : url.slice(slash);
  return path.split('?')[0].split('#')[0];
}

// The pair record at `key` (state/active or state/previous), or null if none.
async function pairAt(bucket, key) {
  const object = await bucket.get(key);
  if (object === null) return null;
  const pair = JSON.parse(await object.text());
  if (!NAME.test(pair.app ?? '') || !NAME.test(pair.publication ?? '')) {
    throw new Error(`${key} names no valid pair`);
  }
  return pair;
}

function etagMatches(header, etag) {
  return header !== null && header.split(',').some(tag => tag.trim().replace(/^W\//, '') === etag || tag.trim() === '*');
}

// serve.rs file_response for the object at `key`; `cache` is its Cache-Control.
async function objectResponse(bucket, key, cache, request) {
  // R2 keys are at most 1,024 bytes, so a longer one names nothing; asking R2
  // for it would fail with 500 rather than answer as a missing file.
  if (new TextEncoder().encode(key).length > 1024) return null;
  const head = request.method === 'HEAD';
  let rangeHeader = request.headers.get('Range');
  const ifRange = request.headers.get('If-Range');
  const conditional = request.headers.get('If-None-Match') !== null;
  // A plain GET, most of the app's requests, takes one read (a Class B
  // operation) instead of a head and a read.
  if (!head && rangeHeader === null && !conditional) {
    const object = await bucket.get(key);
    if (object === null) return null;
    return new Response(object.size ? object.body : null, {
      status: 200,
      headers: {
        'Content-Type': contentType(key),
        'Cache-Control': cache,
        'Accept-Ranges': 'bytes',
        ETag: object.httpEtag,
        'Content-Length': String(object.size),
      },
    });
  }
  const meta = await bucket.head(key);
  if (meta === null) return null;
  if (rangeHeader !== null && ifRange !== null && ifRange !== meta.httpEtag) rangeHeader = null;
  const range = requestedRange(rangeHeader, meta.size);
  if (range === null) {
    return text(416, 'range not satisfiable', { 'Accept-Ranges': 'bytes', 'Content-Range': `bytes */${meta.size}` });
  }
  const headers = new Headers({
    'Content-Type': contentType(key),
    'Cache-Control': cache,
    'Accept-Ranges': 'bytes',
    ETag: meta.httpEtag,
  });
  if (range.full && etagMatches(request.headers.get('If-None-Match'), meta.httpEtag)) {
    return new Response(null, { status: 304, headers });
  }
  const length = range.full ? meta.size : range.length;
  headers.set('Content-Length', String(length));
  if (!range.full) headers.set('Content-Range', `bytes ${range.offset}-${range.offset + length - 1}/${meta.size}`);
  const status = range.full ? 200 : 206;
  if (head || length === 0) return new Response(null, { status, headers });
  const object = await bucket.get(key, range.full ? {} : { range: { offset: range.offset, length } });
  // Releases are immutable, so the object read is the one `head` described.
  if (object === null || object.httpEtag !== meta.httpEtag) throw new Error(`${key} changed while it was served`);
  return new Response(object.body, { status, headers });
}

export async function handle(request, bucket) {
  if (request.method !== 'GET' && request.method !== 'HEAD') {
    return text(405, 'method not allowed', { Allow: 'GET, HEAD' });
  }
  const pair = await pairAt(bucket, ACTIVE_KEY);
  if (pair === null) return text(503, 'no active release');
  const app = `app/${pair.app}/dist/`;
  const index = app + 'index.html';
  const url = rawPath(request.url);
  if (url.startsWith('/data/')) {
    const key = keyUnder(`pub/${pair.publication}/report/`, url.slice('/data/'.length));
    const response = key && (await objectResponse(bucket, key, 'no-cache', request));
    return response || text(404, 'not found');
  }
  const key = url === '/' ? index : keyUnder(app, url);
  if (key !== null) {
    // Hashed build output is immutable; a fallback for a missing asset is not.
    const asset = url.startsWith('/assets/') && key.startsWith(app + 'assets/');
    const response = await objectResponse(bucket, key, asset ? IMMUTABLE : 'no-cache', request);
    if (response) return response;
    // A tab loaded before a switch to a new app label still asks for the old
    // app's chunks. Hashed names can't collide, so serve them from the
    // previous app rather than as index.html, which fails as a module. This
    // is best effort: a state/previous that doesn't parse leaves the miss to
    // index.html, as before, instead of failing it.
    if (asset) {
      let previous = null;
      try {
        previous = await pairAt(bucket, PREVIOUS_KEY);
      } catch (error) {
        console.error(String(error));
      }
      if (previous !== null && previous.app !== pair.app) {
        const old = `app/${previous.app}/dist/` + key.slice(app.length);
        const response = await objectResponse(bucket, old, IMMUTABLE, request);
        if (response) return response;
      }
    }
  }
  const response = await objectResponse(bucket, index, 'no-cache', request);
  if (response) return response;
  throw new Error(`${index} is missing`);
}

export default {
  async fetch(request, env) {
    let response;
    try {
      response = await handle(request, env.SITE);
    } catch (error) {
      console.error(String(error));
      response = text(500, 'internal server error');
    }
    // The runtime drops a HEAD response's body and with it the length of a
    // text answer (404, 405, 416); grm serve keeps it, so set it here.
    if (request.method === 'HEAD' && response.body !== null) {
      const length = (await response.arrayBuffer()).byteLength;
      const headers = new Headers(response.headers);
      headers.set('Content-Length', String(length));
      response = new Response(null, { status: response.status, headers });
    }
    return response;
  },
};
