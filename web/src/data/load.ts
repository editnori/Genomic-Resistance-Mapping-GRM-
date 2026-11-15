// One fetch cache for the report files. Each path is fetched once, parsed at
// the edge, and shared by every view that asks for it.

import { useCallback, useMemo, useSyncExternalStore } from "react";
import type { Corpus, DatasetReport, Manifest, Research, SpeciesReport } from "./schema";
import { ShapeError, parseCorpus, parseDataset, parseManifest, parseResearch, parseSpecies } from "./validate";
import { parseIntervals, type Intervals } from "./intervals";
import { parseBaselines, type Baselines } from "./baselines";

export type LoadError =
  | { kind: "missing"; path: string }
  | { kind: "http"; path: string; status: number }
  | { kind: "network"; path: string; message: string }
  | { kind: "parse"; path: string; message: string }
  | { kind: "schema"; path: string; message: string; found: unknown }
  | { kind: "shape"; path: string; message: string };

export type Resource<T> =
  | { status: "loading" }
  | { status: "ready"; data: T }
  | { status: "error"; error: LoadError };

interface Entry {
  state: Resource<unknown>;
  listeners: Set<() => void>;
  started: boolean;
}

const LOADING: Resource<never> = { status: "loading" };
const cache = new Map<string, Entry>();

const entryFor = (path: string): Entry => {
  let entry = cache.get(path);
  if (entry === undefined) {
    entry = { state: LOADING, listeners: new Set(), started: false };
    cache.set(path, entry);
  }
  return entry;
};

// ---- report digests ------------------------------------------------------------------

interface DigestEntry {
  /** Undefined until the report has been read; null when no digest exists. */
  value: string | null | undefined;
  listeners: Set<() => void>;
}

const digests = new Map<string, DigestEntry>();

const digestEntry = (path: string): DigestEntry => {
  let entry = digests.get(path);
  if (entry === undefined) {
    entry = { value: undefined, listeners: new Set() };
    digests.set(path, entry);
  }
  return entry;
};

const settleDigest = (path: string, value: string | null): void => {
  const entry = digestEntry(path);
  if (entry.value !== undefined) return;
  entry.value = value;
  entry.listeners.forEach((l) => l());
};

/** SHA-256 of bytes as lowercase hex, or null where the page has no Web
 *  Crypto (a page served over plain HTTP from another host than localhost)
 *  or the digest fails. */
export async function sha256Hex(bytes: ArrayBuffer): Promise<string | null> {
  const subtle = globalThis.crypto?.subtle;
  if (subtle === undefined) return null;
  try {
    const digest = await subtle.digest("SHA-256", bytes);
    return Array.from(new Uint8Array(digest), (b) => b.toString(16).padStart(2, "0")).join("");
  } catch {
    return null;
  }
}

/** A dataset report's bytes are hashed as they are read, for the optional
 *  files pinned to them (baselines/<id>.json). */
const isDatasetPath = (path: string): boolean => path.startsWith("/data/datasets/");

/** Raw JSON for a /data path, or the error that stopped it. A dataset
 *  report's SHA-256 is computed from the same bytes, beside the parse, and
 *  settles to null when the report cannot be read. */
async function fetchRaw(path: string): Promise<{ ok: true; raw: unknown } | { ok: false; error: LoadError }> {
  const got = await fetchBytes(path);
  if (isDatasetPath(path)) {
    if (got.ok) void sha256Hex(got.bytes).then((d) => settleDigest(path, d));
    else settleDigest(path, null);
  }
  if (!got.ok) return got;
  try {
    // What response.json() does: UTF-8 without a byte-order mark, then JSON.
    return { ok: true, raw: JSON.parse(new TextDecoder("utf-8").decode(got.bytes)) as unknown };
  } catch (e) {
    return { ok: false, error: { kind: "parse", path, message: e instanceof Error ? e.message : String(e) } };
  }
}

async function fetchBytes(path: string): Promise<{ ok: true; bytes: ArrayBuffer } | { ok: false; error: LoadError }> {
  let response: Response;
  try {
    response = await fetch(path, { headers: { Accept: "application/json" } });
  } catch (e) {
    return { ok: false, error: { kind: "network", path, message: e instanceof Error ? e.message : String(e) } };
  }
  if (response.status === 404) return { ok: false, error: { kind: "missing", path } };
  if (!response.ok) return { ok: false, error: { kind: "http", path, status: response.status } };
  const type = response.headers.get("Content-Type") ?? "";
  // A server that falls back to index.html for unknown paths answers 200
  // with HTML; that is a missing file, not a parse failure.
  if (type.includes("text/html")) return { ok: false, error: { kind: "missing", path } };
  try {
    return { ok: true, bytes: await response.arrayBuffer() };
  } catch (e) {
    return { ok: false, error: { kind: "network", path, message: e instanceof Error ? e.message : String(e) } };
  }
}

async function load<T>(path: string, parse: (raw: unknown) => T): Promise<Resource<T>> {
  const got = await fetchRaw(path);
  if (!got.ok) return { status: "error", error: got.error };
  try {
    return { status: "ready", data: parse(got.raw) };
  } catch (e) {
    if (e instanceof ShapeError) {
      return e.foundSchema !== undefined
        ? { status: "error", error: { kind: "schema", path, message: e.message, found: e.foundSchema } }
        : { status: "error", error: { kind: "shape", path, message: e.message } };
    }
    return { status: "error", error: { kind: "shape", path, message: e instanceof Error ? e.message : String(e) } };
  }
}

function start<T>(path: string, parse: (raw: unknown) => T): void {
  const entry = entryFor(path);
  if (entry.started) return;
  entry.started = true;
  load(path, parse)
    .catch((e: unknown): Resource<T> => ({
      status: "error",
      error: { kind: "network", path, message: e instanceof Error ? e.message : String(e) },
    }))
    .then((state) => {
      entry.state = state;
      entry.listeners.forEach((l) => l());
    });
}

/** The cached state of one report file; starts the fetch on first use. A
 *  null path (nothing selected yet) reads as loading without fetching. */
export function useResource<T>(path: string | null, parse: (raw: unknown) => T): Resource<T> {
  const subscribe = useCallback(
    (onChange: () => void) => {
      if (path === null) return () => undefined;
      const entry = entryFor(path);
      entry.listeners.add(onChange);
      start(path, parse);
      return () => entry.listeners.delete(onChange);
    },
    [path, parse],
  );
  const snapshot = useCallback(() => (path === null ? LOADING : entryFor(path).state), [path]);
  return useSyncExternalStore(subscribe, snapshot, snapshot) as Resource<T>;
}

export const PATHS = {
  manifest: "/data/index.json",
  corpus: "/data/corpus.json",
  species: (slug: string) => `/data/species/${encodeURIComponent(slug)}.json`,
  dataset: (id: string) => `/data/datasets/${encodeURIComponent(id)}.json`,
  research: "/data/research.json",
  intervals: (id: string) => `/data/intervals/${encodeURIComponent(id)}.json`,
  baselines: (id: string) => `/data/baselines/${encodeURIComponent(id)}.json`,
} as const;

export const useManifest = (): Resource<Manifest> => useResource(PATHS.manifest, parseManifest);
export const useCorpus = (): Resource<Corpus> => useResource(PATHS.corpus, parseCorpus);
export const useSpecies = (slug: string | null): Resource<SpeciesReport> =>
  useResource(slug === null ? null : PATHS.species(slug), parseSpecies);
export const useDataset = (id: string | null): Resource<DatasetReport> =>
  useResource(id === null ? null : PATHS.dataset(id), parseDataset);
/** Optional: a missing or invalid research.json is an error state of this
 *  resource only, which researchFor reads as no research. */
export const useResearch = (): Resource<Research> => useResource(PATHS.research, parseResearch);
/** Optional, like research.json: a missing or invalid intervals file reads
 *  as no intervals (intervalsFor). */
export const useIntervals = (id: string): Resource<Intervals> => useResource(PATHS.intervals(id), parseIntervals);
/** Optional too: a missing or invalid baselines file reads as no baselines
 *  (baselinesFor). */
export const useBaselines = (id: string): Resource<Baselines> => useResource(PATHS.baselines(id), parseBaselines);

/** The SHA-256 of dataset `id`'s report as this page read it: undefined
 *  until it has been read and hashed, null when it could not be. */
export function useReportDigest(id: string): string | null | undefined {
  const path = PATHS.dataset(id);
  const subscribe = useCallback(
    (onChange: () => void) => {
      const entry = digestEntry(path);
      entry.listeners.add(onChange);
      return () => entry.listeners.delete(onChange);
    },
    [path],
  );
  const snapshot = useCallback(() => digestEntry(path).value, [path]);
  return useSyncExternalStore(subscribe, snapshot, snapshot);
}

/** Every baselines file of the datasets with an MIC model, by id, read in
 *  one pass so the overview and the table draw once. A file that is missing
 *  or invalid is left out, and that dataset shows no constant. */
async function loadMicBaselines(ids: readonly string[]): Promise<Resource<ReadonlyMap<string, Baselines>>> {
  const files = await Promise.all(
    ids.map(async (id) => {
      const got = await fetchRaw(PATHS.baselines(id));
      if (!got.ok) return null;
      try {
        return parseBaselines(got.raw);
      } catch {
        return null;
      }
    }),
  );
  const out = new Map<string, Baselines>();
  files.forEach((f, i) => {
    if (f !== null) out.set(ids[i]!, f);
  });
  return { status: "ready", data: out };
}

/** The cache key of one manifest's MIC datasets' files: one entry per set
 *  of ids, so another manifest reads its own. */
const micBaselinesKey = (ids: readonly string[]): string => `baselines of ${ids.join(" ")}`;

export function useMicBaselines(manifest: Manifest): Resource<ReadonlyMap<string, Baselines>> {
  const ids = useMemo(() => manifest.datasets.filter((d) => d.mic_essential_agreement !== null).map((d) => d.id), [manifest]);
  const key = micBaselinesKey(ids);
  const subscribe = useCallback(
    (onChange: () => void) => {
      const entry = entryFor(key);
      entry.listeners.add(onChange);
      if (!entry.started) {
        entry.started = true;
        void loadMicBaselines(ids).then((state) => {
          entry.state = state;
          entry.listeners.forEach((l) => l());
        });
      }
      return () => entry.listeners.delete(onChange);
    },
    [key, ids],
  );
  const snapshot = useCallback(() => entryFor(key).state as Resource<ReadonlyMap<string, Baselines>>, [key]);
  return useSyncExternalStore(subscribe, snapshot, snapshot);
}

/** A load failure in words a visitor can act on. Where the cause is a
 *  report that was never built or was built by another version, the last
 *  sentence tells someone who runs GRM what to do. */
export const loadErrorText = (error: LoadError): string => {
  switch (error.kind) {
    case "missing":
      return `${error.path} was not found, so this part of the report is not published here. If you run GRM yourself, \`grm report\` writes this file.`;
    case "http":
      return `${error.path} answered HTTP ${error.status}.`;
    case "network":
      return `${error.path} could not be fetched: ${error.message}.`;
    case "parse":
      return `${error.path} is not valid JSON: ${error.message}.`;
    case "schema":
      return `${error.message}. The page and the report come from different versions of GRM. If you run GRM yourself, rebuild the report with the matching \`grm\` version.`;
    case "shape":
      return `${error.path} does not match the dashboard contract: ${error.message}.`;
  }
};
