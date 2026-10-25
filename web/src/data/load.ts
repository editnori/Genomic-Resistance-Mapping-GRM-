// One fetch cache for the report files. Each path is fetched once, parsed at
// the edge, and shared by every view that asks for it.

import { useCallback, useSyncExternalStore } from "react";
import type { Corpus, DatasetReport, Manifest, Research, SpeciesReport } from "./schema";
import { ShapeError, parseCorpus, parseDataset, parseManifest, parseResearch, parseSpecies } from "./validate";

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

/** Raw JSON for a /data path, or the error that stopped it. */
async function fetchRaw(path: string): Promise<{ ok: true; raw: unknown } | { ok: false; error: LoadError }> {
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
    return { ok: true, raw: await response.json() };
  } catch (e) {
    return { ok: false, error: { kind: "parse", path, message: e instanceof Error ? e.message : String(e) } };
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

/** One literal sentence for a load failure. */
export const loadErrorText = (error: LoadError): string => {
  switch (error.kind) {
    case "missing":
      return `${error.path} does not exist yet. Run \`grm report\` once the pipeline has written this file.`;
    case "http":
      return `${error.path} answered HTTP ${error.status}.`;
    case "network":
      return `${error.path} could not be fetched: ${error.message}.`;
    case "parse":
      return `${error.path} is not valid JSON: ${error.message}.`;
    case "schema":
      return `${error.message}. Rebuild the report with the matching \`grm\` version.`;
    case "shape":
      return `${error.path} does not match the dashboard contract: ${error.message}.`;
  }
};
