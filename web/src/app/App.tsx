// Dark bezel, one rounded paper canvas, the top bar, then one mounted view
// that owns its own scroll. The dataset inspector and Present are overlays
// above whichever view is mounted.

import { lazy, Suspense, useEffect, useRef, useState } from "react";
import { flushSync } from "react-dom";
import { loadErrorText, useManifest } from "../data/load";
import type { Manifest } from "../data/schema";
import { PanelState, VIEW_SHELL } from "../ui/primitives";
import { recordLaunch } from "../ui/launch";
import { TopBar } from "./TopBar";
import { findFocusToken, isAppPath, lastNavKind, navigate, useLocationKey, useRouter, viewScrollOwner, VIEWS, type FocusToken, type Navigator, type ViewId } from "./router";
import { OverviewView } from "../views/OverviewView";
import { DatasetsView } from "../views/DatasetsView";
import { SpeciesView } from "../views/SpeciesView";
import { MethodsView } from "../views/MethodsView";
import { DatasetInspector } from "../views/inspector/DatasetInspector";

const PresentMode = lazy(() => import("../present/PresentMode").then((m) => ({ default: m.PresentMode })));

const ViewBody = ({ manifest, nav }: { manifest: Manifest; nav: Navigator }) => {
  switch (nav.route.view) {
    case "overview":
      return <OverviewView manifest={manifest} nav={nav} />;
    case "datasets":
      return <DatasetsView manifest={manifest} nav={nav} />;
    case "species":
      return <SpeciesView manifest={manifest} nav={nav} />;
    case "methods":
      return <MethodsView manifest={manifest} />;
  }
};

const reducedMotion = (): boolean => window.matchMedia("(prefers-reduced-motion: reduce)").matches;

type ViewTransitionDocument = Document & { startViewTransition?: (update: () => void) => unknown };

/** Switches view inside a view transition where the browser has one: the
 *  old view fades under the new, which steps in through its own class. */
const switchView = (update: () => void): void => {
  const doc = document as ViewTransitionDocument;
  if (typeof doc.startViewTransition !== "function" || reducedMotion()) {
    update();
    return;
  }
  doc.startViewTransition(() => flushSync(update));
};

/** The step a newly mounted view takes: from the side the reader travelled
 *  toward in the tab order; none on the first view of a visit. */
const useViewStep = (view: ViewId): string => {
  const index = VIEWS.findIndex((v) => v.id === view);
  const [shown, setShown] = useState({ index, step: "" });
  if (shown.index === index) return shown.step;
  // A new view: remember which way the reader went until the next switch.
  const step = index > shown.index ? "view-step-forward" : "view-step-back";
  setShown({ index, step });
  return step;
};

/** Focus that lands on an element the reader cannot tab to (a heading, a
 *  section) without scrolling it again. */
const focusTarget = (el: HTMLElement): void => {
  if (!el.matches("a[href], button, input, select, textarea, [tabindex]")) {
    el.setAttribute("tabindex", "-1");
    el.setAttribute("data-focus-target", "");
  }
  el.focus({ preventScroll: true });
};

const viewHeading = (): HTMLElement | null => viewScrollOwner()?.querySelector<HTMLElement>("h1") ?? null;

/** A modal (the inspector, Present, the filter sheet) owns focus. */
const modalOpen = (): boolean => document.querySelector("dialog[open]") !== null;

const focusLost = (): boolean => {
  const a = document.activeElement;
  return a === null || a === document.body || a === document.documentElement || !a.isConnected;
};

/** Reveals a hash's section in the view's scroll owner and moves focus to
 *  its heading. A section drawn later (data still loading) is awaited for a
 *  few seconds; a reader who scrolls or types first keeps their place. */
const revealHash = (hash: string): (() => void) => {
  const id = decodeURIComponent(hash.replace(/^#/, ""));
  let done = false;
  const reveal = (): boolean => {
    const el = document.getElementById(id);
    if (el === null) return false;
    el.scrollIntoView({ block: "start" });
    const heading = el.matches("h1, h2, h3") ? el : el.querySelector<HTMLElement>("h1, h2, h3");
    focusTarget(heading ?? el);
    return true;
  };
  if (reveal()) return () => {};
  const stop = () => {
    done = true;
    observer.disconnect();
    window.clearTimeout(timer);
    window.removeEventListener("wheel", stop, true);
    window.removeEventListener("keydown", stop, true);
    window.removeEventListener("pointerdown", stop, true);
  };
  const observer = new MutationObserver(() => {
    if (!done && reveal()) stop();
  });
  observer.observe(document.body, { childList: true, subtree: true });
  const timer = window.setTimeout(stop, 6000);
  window.addEventListener("wheel", stop, true);
  window.addEventListener("keydown", stop, true);
  window.addEventListener("pointerdown", stop, true);
  return stop;
};

/** Back to an entry: its scroll offset once the view is tall enough, then
 *  the control that left it, else the view's title. */
const restoreEntry = (scroll: number | undefined, focus: FocusToken | undefined): (() => void) => {
  let frame = 0;
  let tries = 0;
  const step = () => {
    const owner = viewScrollOwner();
    if (owner !== null && scroll !== undefined && scroll > 0) {
      owner.scrollTop = scroll;
      // A view whose data is still arriving is too short for the offset yet.
      if (Math.abs(owner.scrollTop - scroll) > 1 && tries < 120) {
        tries += 1;
        frame = window.requestAnimationFrame(step);
        return;
      }
    }
    if (modalOpen()) return;
    const el = focus === undefined ? null : findFocusToken(focus);
    if (el !== null) el.focus({ preventScroll: true });
    else if (focusLost()) {
      const h1 = viewHeading();
      if (h1 !== null) focusTarget(h1);
    }
  };
  step();
  return () => window.cancelAnimationFrame(frame);
};

export const App = () => {
  const nav = useRouter();
  const at = useLocationKey();
  const manifest = useManifest();
  const [presenting, setPresenting] = useState(false);
  const step = useViewStep(nav.route.view);
  const ready = manifest.status === "ready";

  // The inspector opens from the row or mark that launched it.
  useEffect(() => {
    const launcher = (target: EventTarget | null) =>
      target instanceof Element ? (target.closest("tr, [data-launch], li, [role='option'], button, a") ?? target) : null;
    const onClick = (e: MouseEvent) => recordLaunch(launcher(e.target));
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Enter") recordLaunch(launcher(e.target));
    };
    document.addEventListener("click", onClick, true);
    document.addEventListener("keydown", onKey, true);
    return () => {
      document.removeEventListener("click", onClick, true);
      document.removeEventListener("keydown", onKey, true);
    };
  }, []);

  // A plain click on an in-app link navigates in place, as the views' own
  // controls do; a modified click, a download or another target is the
  // browser's.
  useEffect(() => {
    const onClick = (e: MouseEvent) => {
      if (e.defaultPrevented || e.button !== 0 || e.metaKey || e.ctrlKey || e.shiftKey || e.altKey) return;
      const link = e.target instanceof Element ? e.target.closest<HTMLAnchorElement>("a[href]") : null;
      if (link === null || (link.target !== "" && link.target !== "_self") || link.hasAttribute("download")) return;
      const url = new URL(link.href, window.location.href);
      if (url.origin !== window.location.origin || !isAppPath(url.pathname)) return;
      e.preventDefault();
      const href = `${url.pathname}${url.search}${url.hash}`;
      if (url.pathname === window.location.pathname) navigate(href);
      else switchView(() => navigate(href));
    };
    document.addEventListener("click", onClick);
    return () => document.removeEventListener("click", onClick);
  }, []);

  // After each move: a hash reveals its section; Back restores the entry's
  // scroll and the control that left it; a new view takes focus at its
  // title; a closed inspector hands focus back to its dataset's row.
  const was = useRef<{ view: ViewId; path: string; hash: string; dataset: string | null } | null>(null);
  useEffect(() => {
    if (!ready) return undefined;
    const now = { view: nav.route.view, path: window.location.pathname, hash: window.location.hash, dataset: nav.route.dataset };
    const before = was.current;
    was.current = now;
    const kind = lastNavKind();
    const state: unknown = window.history.state;
    const kept = state !== null && typeof state === "object" ? (state as { scroll?: number; focus?: FocusToken }) : {};

    if (before === null) return now.hash !== "" ? revealHash(now.hash) : undefined;

    if (before.dataset !== null && now.dataset === null && before.path === now.path) {
      if (kind === "pop" && kept.focus !== undefined) return restoreEntry(undefined, kept.focus);
      if (focusLost()) {
        const row = findFocusToken({ row: before.dataset }, viewScrollOwner() ?? document);
        const target = row ?? viewHeading();
        if (target !== null) {
          if (row !== null) row.scrollIntoView({ block: "nearest" });
          focusTarget(target);
        }
      }
      return undefined;
    }
    if (now.dataset !== null && before.dataset === null && before.path === now.path) return undefined;

    if (kind === "pop") return restoreEntry(kept.scroll, kept.focus);
    if (now.hash !== "" && (now.hash !== before.hash || now.path !== before.path || kind === "replace")) return revealHash(now.hash);
    if (now.view !== before.view && !modalOpen()) {
      const h1 = viewHeading();
      if (h1 !== null) focusTarget(h1);
    }
    return undefined;
  }, [at, ready, nav.route.view, nav.route.dataset]);

  return (
    <div className="flex h-full min-h-0 flex-col bg-bezel p-1 md:p-2">
      <div className="flex min-h-0 flex-1 flex-col overflow-hidden rounded-[14px] bg-paper text-ink shadow-[var(--surface-shadow)] md:rounded-[18px]">
        <TopBar
          active={nav.route.view}
          loading={manifest.status === "loading"}
          onNavigate={(view) => switchView(() => nav.go({ view, species: null }))}
          onPresent={() => setPresenting(true)}
          presentDisabled={manifest.status !== "ready"}
        />
        <main className="relative min-h-0 flex-1 border-t border-hair">
          {/* relative: a visually hidden label anywhere in the view is
              positioned inside this scroll owner, never against main, so it
              cannot make the document itself scroll. */}
          <div key={nav.route.view} data-view-scroll className={`relative h-full overflow-y-auto overflow-x-hidden overscroll-contain ${step}`}>
            {manifest.status === "ready" ? (
              <ViewBody manifest={manifest.data} nav={nav} />
            ) : (
              <div className={`${VIEW_SHELL} py-10`}>
                {manifest.status === "loading" ? (
                  <PanelState kind="loading" title="Loading the report manifest" />
                ) : (
                  <PanelState kind="error" title="The report manifest could not be read" message={loadErrorText(manifest.error)} />
                )}
              </div>
            )}
          </div>
        </main>
      </div>
      {manifest.status === "ready" && nav.route.dataset !== null ? (
        <DatasetInspector key={nav.route.dataset} id={nav.route.dataset} manifest={manifest.data} nav={nav} />
      ) : null}
      {presenting && manifest.status === "ready" ? (
        <Suspense fallback={null}>
          <PresentMode manifest={manifest.data} onClose={() => setPresenting(false)} />
        </Suspense>
      ) : null}
    </div>
  );
};
