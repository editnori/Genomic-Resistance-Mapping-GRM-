// URL routing without a library. Paths name the view (and the species);
// `?dataset=<id>` opens the dataset inspector over whichever view is shown,
// so a dataset link works from any view and Back closes the inspector.

import { useCallback, useSyncExternalStore } from "react";

export const VIEWS = [
  { id: "overview", label: "Overview", path: "/" },
  { id: "datasets", label: "Datasets", path: "/datasets" },
  { id: "species", label: "Species", path: "/species" },
  { id: "methods", label: "Data & Methods", path: "/methods" },
] as const;

/** Paths that moved, so old links keep working: the Molecular tests view is
 *  now the panel block's evidence on Overview. */
const MOVED: Record<string, { path: string; hash: string }> = { "/tests": { path: "/", hash: "#panels" } };

/** Where a moved path now lives, keeping its query (an open dataset); null
 *  for a path that did not move. */
export const movedHref = (pathname: string, search: string): string | null => {
  const first = pathname.split("/").filter((p) => p !== "")[0];
  const to = first === undefined ? undefined : MOVED[`/${first}`];
  return to === undefined ? null : `${to.path}${search}${to.hash}`;
};

/** Replaces a moved path before the first render, so Back skips it. */
export const redirectMovedPath = (): void => {
  const href = movedHref(window.location.pathname, window.location.search);
  if (href !== null) window.history.replaceState(null, "", href);
};

export type ViewId = (typeof VIEWS)[number]["id"];

export interface Route {
  view: ViewId;
  /** /species/<slug> */
  species: string | null;
  /** ?dataset=<id> */
  dataset: string | null;
}

export const parseRoute = (pathname: string, search: string): Route => {
  const parts = pathname.split("/").filter((p) => p !== "").map(decodeURIComponent);
  const view = (VIEWS.find((v) => v.path === `/${parts[0] ?? ""}`)?.id ?? "overview") as ViewId;
  return {
    view,
    species: view === "species" ? (parts[1] ?? null) : null,
    dataset: new URLSearchParams(search).get("dataset"),
  };
};

export const routeHref = (r: Route): string => {
  const base = VIEWS.find((v) => v.id === r.view)!.path;
  const path = r.view === "species" && r.species !== null ? `${base}/${encodeURIComponent(r.species)}` : base;
  return r.dataset === null ? path : `${path}?dataset=${encodeURIComponent(r.dataset)}`;
};

const listeners = new Set<() => void>();
const notify = () => listeners.forEach((l) => l());
const subscribe = (l: () => void) => {
  listeners.add(l);
  window.addEventListener("popstate", l);
  return () => {
    listeners.delete(l);
    window.removeEventListener("popstate", l);
  };
};
const snapshot = () => `${window.location.pathname}${window.location.search}`;

interface HistoryMark {
  /** Set on entries pushed by opening the inspector, so closing can go Back. */
  inspector?: boolean;
}

export interface Navigator {
  route: Route;
  go(next: Partial<Route>): void;
  openDataset(id: string): void;
  closeDataset(): void;
}

export const useRouter = (): Navigator => {
  const href = useSyncExternalStore(subscribe, snapshot, snapshot);
  const url = new URL(href, "http://local");
  const route = parseRoute(url.pathname, url.search);

  const go = useCallback((next: Partial<Route>) => {
    const cur = parseRoute(window.location.pathname, window.location.search);
    const target = { ...cur, dataset: null, ...next };
    window.history.pushState({} satisfies HistoryMark, "", routeHref(target));
    notify();
  }, []);

  const openDataset = useCallback((id: string) => {
    const cur = parseRoute(window.location.pathname, window.location.search);
    const mark: HistoryMark = { inspector: true };
    if (cur.dataset !== null) window.history.replaceState(mark, "", routeHref({ ...cur, dataset: id }));
    else window.history.pushState(mark, "", routeHref({ ...cur, dataset: id }));
    notify();
  }, []);

  const closeDataset = useCallback(() => {
    const mark = window.history.state as HistoryMark | null;
    if (mark?.inspector === true) {
      window.history.back();
      return;
    }
    const cur = parseRoute(window.location.pathname, window.location.search);
    window.history.replaceState({}, "", routeHref({ ...cur, dataset: null }));
    notify();
  }, []);

  return { route, go, openDataset, closeDataset };
};
