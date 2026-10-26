// Dark bezel, one rounded paper canvas, the top bar, then one mounted view
// that owns its own scroll. The dataset inspector and Present are overlays
// above whichever view is mounted.

import { lazy, Suspense, useEffect, useState } from "react";
import { flushSync } from "react-dom";
import { loadErrorText, useManifest } from "../data/load";
import type { Manifest } from "../data/schema";
import { PanelState, VIEW_SHELL } from "../ui/primitives";
import { recordLaunch } from "../ui/launch";
import { TopBar } from "./TopBar";
import { useRouter, VIEWS, type Navigator, type ViewId } from "./router";
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

export const App = () => {
  const nav = useRouter();
  const manifest = useManifest();
  const [presenting, setPresenting] = useState(false);
  const step = useViewStep(nav.route.view);

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
