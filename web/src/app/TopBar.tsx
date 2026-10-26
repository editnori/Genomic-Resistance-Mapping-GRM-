// The product's one navigation row: the mark and wordmark (a link to
// Overview), one rounded pill group with one neutral active pill, then
// Present and the theme toggle. When the pills cannot fit (measured, not
// guessed), the compact view selector takes their place; type never shrinks
// to keep them. The brand shares the views' 24px left edge.

import { useLayoutEffect, useRef, useState, type MouseEvent } from "react";
import { Play } from "lucide-react";
import { ThemeToggle } from "../theme/ThemeProvider";
import { Dropdown } from "../ui/Dropdown";
import { FOCUS_RING } from "../ui/primitives";
import { GrmMark, useMarkMotion } from "./Mark";
import { VIEWS, type ViewId } from "./router";

const FIT_ALLOWANCE_PX = 96;

/** A plain click navigates in place; a modified click opens a new tab. */
const plainClick = (e: MouseEvent): boolean => e.button === 0 && !e.metaKey && !e.ctrlKey && !e.shiftKey && !e.altKey;

export const TopBar = ({
  active,
  loading,
  onNavigate,
  onPresent,
  presentDisabled,
}: {
  active: ViewId;
  /** The report is still loading: the mark runs its loop. */
  loading: boolean;
  onNavigate(view: ViewId): void;
  onPresent(): void;
  presentDisabled: boolean;
}) => {
  const headerRef = useRef<HTMLElement | null>(null);
  const brandRef = useRef<HTMLAnchorElement | null>(null);
  const rightRef = useRef<HTMLDivElement | null>(null);
  const measureRef = useRef<HTMLDivElement | null>(null);
  const [fits, setFits] = useState(false);
  const mark = useMarkMotion(loading);

  useLayoutEffect(() => {
    const update = () => {
      const [h, b, r, m] = [headerRef.current, brandRef.current, rightRef.current, measureRef.current];
      if (h === null || b === null || r === null || m === null) return;
      // The pills sit centered, so they must clear the wider side twice.
      setFits(m.offsetWidth + Math.max(b.offsetWidth, r.offsetWidth) * 2 + FIT_ALLOWANCE_PX <= h.clientWidth);
    };
    update();
    const observer = new ResizeObserver(update);
    [headerRef, measureRef, rightRef, brandRef].forEach((ref) => ref.current !== null && observer.observe(ref.current));
    return () => observer.disconnect();
  }, []);

  return (
    <header
      ref={headerRef}
      className={`relative grid h-14 shrink-0 items-center bg-paper px-4 sm:px-6 ${
        fits ? "grid-cols-[minmax(0,1fr)_auto_minmax(0,1fr)] gap-4" : "grid-cols-[auto_minmax(0,1fr)_auto] gap-2"
      }`}
    >
      <a
        ref={brandRef}
        href="/"
        onClick={(e) => {
          if (!plainClick(e)) return;
          e.preventDefault();
          onNavigate("overview");
        }}
        onPointerEnter={mark.replay}
        onFocus={mark.replay}
        className={`-mx-1 flex min-w-0 items-center gap-2.5 justify-self-start rounded-md px-1 py-1 text-ink ${FOCUS_RING}`}
      >
        <GrmMark size={24} motion={mark.motion} />
        <span className="whitespace-nowrap text-body font-semibold uppercase tracking-widest">
          <span className="hidden sm:inline">Genomic Resistance Mapping</span>
          <span className="sm:hidden">GRM</span>
        </span>
      </a>

      {/* Invisible measuring copy at the active weight, so selecting a pill
          never flips the fit. */}
      <div ref={measureRef} aria-hidden="true" className="pointer-events-none invisible absolute left-0 top-0 flex items-center gap-0.5 whitespace-nowrap rounded-full border border-hair p-1">
        {VIEWS.map((v) => (
          <span key={v.id} className="rounded-full px-3 py-1.5 text-control font-semibold leading-none">
            {v.label}
          </span>
        ))}
      </div>

      {fits ? (
        <nav aria-label="Views" className="justify-self-center">
          <div className="flex items-center gap-0.5 rounded-full border border-hair bg-paper-pure p-1">
            {VIEWS.map((v) => {
              const current = v.id === active;
              return (
                <button
                  key={v.id}
                  type="button"
                  onClick={() => onNavigate(v.id)}
                  aria-current={current ? "page" : undefined}
                  className={`rounded-full px-3 py-1.5 text-control leading-none transition-[background-color,color,transform] duration-fast ease-out active:scale-[0.98] ${FOCUS_RING} ${
                    current ? "bg-surface-control font-semibold text-ink" : "text-muted-text hover:bg-surface-subtle hover:text-ink"
                  }`}
                >
                  {v.label}
                </button>
              );
            })}
          </div>
        </nav>
      ) : (
        <nav aria-label="View" className="min-w-0 max-w-full justify-self-center">
          <Dropdown ariaLabel="View" options={VIEWS.map((v) => ({ id: v.id, label: v.label }))} mode={{ kind: "one", value: active, onSelect: onNavigate }} />
        </nav>
      )}

      <div ref={rightRef} className="flex items-center gap-2 justify-self-end sm:gap-3">
        <button
          type="button"
          onClick={onPresent}
          disabled={presentDisabled}
          aria-label="Present"
          className={`ui-control-target inline-flex h-8 items-center gap-1.5 rounded-full border border-hair bg-paper-pure px-2.5 text-control font-medium leading-none text-ink transition-[color,border-color,transform] duration-fast ease-out active:scale-[0.98] enabled:hover:border-blue enabled:hover:text-blue disabled:cursor-not-allowed disabled:opacity-45 min-[1280px]:px-4 ${FOCUS_RING}`}
        >
          <Play aria-hidden="true" className="h-4 w-4" strokeWidth={1.75} />
          <span className="hidden min-[1280px]:inline">Present</span>
        </button>
        <ThemeToggle />
      </div>
    </header>
  );
};
