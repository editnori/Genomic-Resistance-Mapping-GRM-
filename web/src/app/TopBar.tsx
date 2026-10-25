import { Play } from "lucide-react";
import type { Manifest } from "../data/schema";
import { ThemeToggle } from "../theme/ThemeProvider";
import { FOCUS_RING } from "../ui/primitives";
import { Brand } from "./Brand";
import { ReportDetails } from "./ReportDetails";
import { VIEWS, type ViewId } from "./router";

export const TopBar = ({ active, onNavigate, onPresent, manifest }: {
  active: ViewId;
  onNavigate(view: ViewId): void;
  onPresent(): void;
  manifest?: Manifest;
}) => (
  <header className="grid shrink-0 grid-cols-[1fr_auto] items-center gap-x-4 bg-paper px-4 sm:px-6 min-[900px]:grid-cols-[1fr_auto_1fr]">
    <button type="button" onClick={() => onNavigate("overview")} aria-label="GRM overview" className={`flex min-h-14 w-fit items-center rounded-md ${FOCUS_RING}`}><Brand /></button>
    <nav aria-label="Views" className="order-3 col-span-2 flex min-h-11 items-center gap-1 overflow-x-auto min-[900px]:order-none min-[900px]:col-span-1">
      {VIEWS.map(view => <button key={view.id} type="button" onClick={() => onNavigate(view.id)} aria-current={view.id === active ? "page" : undefined}
        className={`min-h-11 whitespace-nowrap rounded-md px-3 text-control ${FOCUS_RING} ${view.id === active ? "bg-surface-control font-semibold text-ink" : "text-muted-text hover:bg-surface-subtle hover:text-ink"}`}>{view.label}</button>)}
    </nav>
    <div className="flex items-center justify-self-end gap-1">
      {manifest !== undefined ? <ReportDetails manifest={manifest} /> : null}
      <button type="button" onClick={onPresent} disabled={manifest === undefined} aria-label="Present"
        className={`inline-flex min-h-11 min-w-11 items-center justify-center gap-1.5 rounded-md px-3 text-control text-ink hover:bg-surface-control disabled:cursor-not-allowed disabled:opacity-45 ${FOCUS_RING}`}>
        <Play size={16} aria-hidden="true" /><span className="hidden sm:inline">Present</span>
      </button>
      <ThemeToggle />
    </div>
  </header>
);
