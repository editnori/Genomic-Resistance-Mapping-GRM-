// The top bar, then one mounted view that owns its own scroll, on one paper
// surface with no bezel. The dataset inspector and Present are overlays
// above whichever view is mounted.

import { lazy, Suspense, useState } from "react";
import { loadErrorText, useManifest } from "../data/load";
import type { Manifest } from "../data/schema";
import { PanelState, WIDE_SHELL } from "../ui/primitives";
import { TopBar } from "./TopBar";
import { useRouter, type Navigator } from "./router";
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

export const App = () => {
  const nav = useRouter();
  const manifest = useManifest();
  const [presenting, setPresenting] = useState(false);

  return (
    <div className="flex h-full min-h-0 flex-col bg-paper text-ink">
      <TopBar
        active={nav.route.view}
        onNavigate={(view) => nav.go({ view, species: null })}
        onPresent={() => setPresenting(true)}
        manifest={manifest.status === "ready" ? manifest.data : undefined}
      />
      <main className="relative min-h-0 flex-1 border-t border-hair">
        <div key={nav.route.view} data-view-scroll className="view-in h-full overflow-y-auto overflow-x-hidden overscroll-contain">
          {manifest.status === "ready" ? (
            <ViewBody manifest={manifest.data} nav={nav} />
          ) : (
            <div className={`${WIDE_SHELL} py-10`}>
              {manifest.status === "loading" ? (
                <PanelState kind="loading" title="Loading the report manifest" />
              ) : (
                <PanelState kind="error" title="The report manifest could not be read" message={loadErrorText(manifest.error)} />
              )}
            </div>
          )}
        </div>
      </main>
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
