// URL routing without a library. Paths name the view (and the species);
// `?dataset=<id>` opens the dataset inspector over whichever view is shown,
// so a dataset link works from any view and Back closes the inspector. A
// hash names a section; opening and closing the inspector keep it. Each
// history entry also keeps what the reader left there (a view's own state,
// the scroll position, the control in focus), so Back returns them to it.

import { useCallback, useEffect, useRef, useState, useSyncExternalStore, type Dispatch, type SetStateAction } from "react";

export const VIEWS = [
  { id: "overview", label: "Overview", path: "/" },
  { id: "datasets", label: "Datasets", path: "/datasets" },
  { id: "species", label: "Species", path: "/species" },
  { id: "methods", label: "Data & Methods", path: "/methods" },
] as const;

/** The Overview section the status line's audit words link to. */
export const LIMITS_ID = "limits";
export const LIMITS_HREF = `/#${LIMITS_ID}`;

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

/** A hash as the URL writes it: "" or "#name". */
const hashOf = (hash: string | undefined): string => (hash === undefined || hash === "" || hash === "#" ? "" : hash.startsWith("#") ? hash : `#${hash}`);

export const routeHref = (r: Route, hash?: string): string => {
  const base = VIEWS.find((v) => v.id === r.view)!.path;
  const path = r.view === "species" && r.species !== null ? `${base}/${encodeURIComponent(r.species)}` : base;
  return `${r.dataset === null ? path : `${path}?dataset=${encodeURIComponent(r.dataset)}`}${hashOf(hash)}`;
};

/** A path the app draws, as opposed to a data file or an asset. */
export const isAppPath = (pathname: string): boolean => {
  const first = pathname.split("/").filter((p) => p !== "")[0];
  if (first === undefined) return true;
  if (MOVED[`/${first}`] !== undefined) return true;
  return VIEWS.some((v) => v.path === `/${first}`) && (first === "species" || pathname.split("/").filter((p) => p !== "").length === 1);
};

// ---- What a history entry keeps -------------------------------------------

/** The control in focus when the reader left an entry, by what finds it
 *  again after the view is drawn anew: a table row's id, an element id, a
 *  link's href, or a control's name. */
export type FocusToken = { row: string } | { id: string } | { href: string } | { label: string };

interface EntryState {
  /** Set on entries pushed by opening the inspector, so closing can go Back. */
  inspector?: boolean;
  /** Each mounted view's own state, by key. */
  views?: Record<string, unknown>;
  /** The view's scroll offset when the reader left the entry. */
  scroll?: number;
  focus?: FocusToken;
}

const entry = (): EntryState => {
  // A server render has no history: every view starts as it arrives.
  const s: unknown = typeof window === "undefined" ? null : window.history.state;
  return s !== null && typeof s === "object" ? (s as EntryState) : {};
};

/** The element that owns a view's vertical scroll. */
export const viewScrollOwner = (): HTMLElement | null => document.querySelector<HTMLElement>("[data-view-scroll]");

const controlName = (el: Element): string => (el.getAttribute("aria-label") ?? el.textContent ?? "").replace(/\s+/g, " ").trim();

export const focusToken = (el: Element | null): FocusToken | null => {
  if (!(el instanceof HTMLElement) || el === document.body || el === document.documentElement) return null;
  const row = el.closest<HTMLElement>("[data-row-open]");
  if (row?.dataset.rowKey !== undefined) return { row: row.dataset.rowKey };
  if (el.id !== "") return { id: el.id };
  const link = el.closest("a[href]");
  if (link !== null) return { href: link.getAttribute("href")! };
  const name = controlName(el);
  return name === "" ? null : { label: name };
};

export const findFocusToken = (t: FocusToken, root: ParentNode = document): HTMLElement | null => {
  if ("row" in t) return [...root.querySelectorAll<HTMLElement>("[data-row-open]")].find((e) => e.dataset.rowKey === t.row) ?? null;
  if ("id" in t) return document.getElementById(t.id);
  if ("href" in t) return [...root.querySelectorAll<HTMLElement>("a[href]")].find((a) => a.getAttribute("href") === t.href) ?? null;
  return [...root.querySelectorAll<HTMLElement>("button, a[href], [role='button'], summary")].find((b) => controlName(b) === t.label) ?? null;
};

/** Views register a reader of their own state; leaving an entry saves it. */
const keepers = new Map<string, () => unknown>();

/** Writes what the reader leaves behind into the current entry, just
 *  before a new entry is pushed over it. */
const saveEntry = (): void => {
  const views: Record<string, unknown> = { ...entry().views };
  keepers.forEach((read, key) => {
    views[key] = read();
  });
  const focus = focusToken(document.activeElement);
  try {
    window.history.replaceState({ ...entry(), views, scroll: viewScrollOwner()?.scrollTop ?? 0, focus: focus ?? undefined } satisfies EntryState, "", window.location.href);
  } catch {
    // A value history cannot clone, or a throttled history: the entry keeps
    // what it had, and Back still works without the extras.
  }
};

/** The state a view kept in this entry, if any. */
export const keptState = <T,>(key: string): T | undefined => entry().views?.[key] as T | undefined;

/** A view's own state that survives leaving and coming Back: it starts from
 *  what this entry kept, else from `arrive`, and is saved into the entry
 *  when the reader moves on. Values must be structured-cloneable. */
export const useKeptState = <T,>(key: string, arrive: T): [T, Dispatch<SetStateAction<T>>] => {
  const [value, setValue] = useState<T>(() => keptState<T>(key) ?? arrive);
  const latest = useRef(value);
  useEffect(() => {
    latest.current = value;
  }, [value]);
  useEffect(() => {
    const read = () => latest.current;
    keepers.set(key, read);
    return () => {
      if (keepers.get(key) === read) keepers.delete(key);
    };
  }, [key]);
  return [value, setValue];
};

// ---- Location store --------------------------------------------------------

/** How the location last changed: the first load, a new entry, an entry
 *  rewritten in place, or Back/Forward. */
export type NavKind = "load" | "push" | "replace" | "pop";

let lastKind: NavKind = "load";
let serial = 0;
const listeners = new Set<() => void>();

const changed = (kind: NavKind): void => {
  lastKind = kind;
  serial += 1;
  listeners.forEach((l) => l());
};

const onPop = (): void => changed("pop");

const subscribe = (l: () => void) => {
  if (listeners.size === 0) window.addEventListener("popstate", onPop);
  listeners.add(l);
  return () => {
    listeners.delete(l);
    if (listeners.size === 0) window.removeEventListener("popstate", onPop);
  };
};

const snapshot = () => (typeof window === "undefined" ? `${serial}|/` : `${serial}|${window.location.pathname}${window.location.search}${window.location.hash}`);

export const lastNavKind = (): NavKind => lastKind;

/** Changes on every location change, including one that only re-reveals
 *  the same hash. */
export const useLocationKey = (): string => useSyncExternalStore(subscribe, snapshot, snapshot);

/** Follows an in-app href: a new entry for a new place (or the same place
 *  with a new hash). The same href again only re-reveals its section. False
 *  for an href the app does not draw. */
export const navigate = (href: string): boolean => {
  const url = new URL(href, window.location.href);
  if (url.origin !== window.location.origin || !isAppPath(url.pathname)) return false;
  const target = movedHref(url.pathname, url.search) ?? `${url.pathname}${url.search}${url.hash}`;
  if (target === `${window.location.pathname}${window.location.search}${window.location.hash}`) {
    changed("replace");
    return true;
  }
  saveEntry();
  window.history.pushState({} satisfies EntryState, "", target);
  changed("push");
  return true;
};

export interface Navigator {
  route: Route;
  /** A new place; `hash` names a section in it. */
  go(next: Partial<Route> & { hash?: string }): void;
  openDataset(id: string): void;
  closeDataset(): void;
}

export const useRouter = (): Navigator => {
  const key = useLocationKey();
  const url = new URL(key.slice(key.indexOf("|") + 1), "http://local");
  const route = parseRoute(url.pathname, url.search);

  const go = useCallback((next: Partial<Route> & { hash?: string }) => {
    const { hash, ...place } = next;
    const cur = parseRoute(window.location.pathname, window.location.search);
    saveEntry();
    window.history.pushState({} satisfies EntryState, "", routeHref({ ...cur, dataset: null, ...place }, hash));
    changed("push");
  }, []);

  const openDataset = useCallback((id: string) => {
    const cur = parseRoute(window.location.pathname, window.location.search);
    const href = routeHref({ ...cur, dataset: id }, window.location.hash);
    if (cur.dataset !== null) {
      window.history.replaceState({ ...entry(), inspector: true } satisfies EntryState, "", href);
      changed("replace");
      return;
    }
    saveEntry();
    // The inspector's entry carries the view's state too, so a view drawn
    // anew under it (Back from another view) comes back the same.
    window.history.pushState({ inspector: true, views: entry().views } satisfies EntryState, "", href);
    changed("push");
  }, []);

  const closeDataset = useCallback(() => {
    const mark = entry();
    if (mark.inspector === true) {
      window.history.back();
      return;
    }
    const cur = parseRoute(window.location.pathname, window.location.search);
    window.history.replaceState({ ...mark, inspector: undefined } satisfies EntryState, "", routeHref({ ...cur, dataset: null }, window.location.hash));
    changed("replace");
  }, []);

  return { route, go, openDataset, closeDataset };
};
