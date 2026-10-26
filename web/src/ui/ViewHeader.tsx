import type { ReactNode } from "react";

/** Every view opens the same way: the visible title, then one quiet status
 *  line of inline facts. A purpose line shows only where the title alone
 *  would leave a reader guessing, never one that restates the title. */
export const ViewHeader = ({ title, purpose, actions, status }: { title: ReactNode; purpose?: string; actions?: ReactNode; status: ReactNode }) => (
  <header>
    <div className="flex flex-wrap items-start justify-between gap-x-8 gap-y-3">
      <div className="min-w-0">
        <h1 className="text-display font-semibold text-ink">{title}</h1>
        {purpose !== undefined ? <p className="type-body mt-2">{purpose}</p> : null}
      </div>
      {actions !== undefined ? <div className="flex min-w-0 flex-wrap items-center gap-2">{actions}</div> : null}
    </div>
    <div className="mt-3 flex min-w-0 flex-wrap items-center gap-x-2 gap-y-1 text-note text-muted-text">{status}</div>
  </header>
);
