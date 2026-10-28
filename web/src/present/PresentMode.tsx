// Present: one deck on a full canvas. A title slide where the mark runs its
// full choreography, the five primer steps (one real report read from zero),
// the five Overview findings, then the rule showcase.
// Each slide stages its beats (the answer, the figure, the caption) and the
// footer's ticks show them land; advancing mid-reveal completes the slide
// first. Only the deck's own controls move it: Previous and Next, the arrow
// and page keys, and the outline for a direct jump; a click on a slide does
// nothing. A figure too dense for the canvas at its minimum row height
// continues onto further slides that repeat its title, answer, and caption,
// so nothing scrolls and no label shrinks. An off-screen copy of each such
// slide measures how many pages it needs. The header, slide, and footer
// share one column. A fixed figure sits one set gap under the answer at
// every width rather than centred in what is left, so a short slide has no
// band of empty space above its figure. No title or outline label is cut
// off: on a phone they wrap, and the footer puts the outline on its own row
// above Previous and Next.

import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState, type MutableRefObject, type ReactNode } from "react";
import type { Manifest } from "../data/schema";
import { dateText, pct, speciesShort } from "../data/format";
import { manifestAudit } from "../data/derive/evaluation";
import { GrmMark } from "../app/Mark";
import { ModalDialog } from "../ui/ModalDialog";
import { Dropdown } from "../ui/Dropdown";
import { ACTION_BUTTON, type Fit, type SlidePart } from "../ui/primitives";
import { FINDINGS, type PagedFit } from "../views/findings";
import { Reveal, useSlideTimeline, type SlideTimeline } from "./motion";
import { RuleFigure, ruleDataset, ruleLineageScore } from "./figures";
import { primerSteps } from "../views/primer";

const COLUMN = "mx-auto w-full max-w-[72rem] px-4 sm:px-8";
/** A fixed figure's box, live and measuring alike: top-aligned, one gap
 *  under the answer. */
const FIXED_BOX = "flex min-h-0 flex-1 flex-col justify-start pb-4 pt-6";
const KEY_OWNERS = 'input, select, textarea, [role="combobox"], [role="listbox"], [role="option"], [contenteditable]:not([contenteditable="false"])';

const owns = (target: EventTarget | null, selector: string): boolean => target instanceof Element && target.closest(selector) !== null;

/** Which page of a slide to draw; `pages` is undefined while measuring. */
interface Paging {
  page: number;
  pages: number | undefined;
}

interface FigureSlide {
  kind: "figure";
  id: string;
  title: string;
  result?: string;
  caption: string;
  /** A figure that fills the height the text leaves. */
  fill?: (fit: PagedFit) => ReactNode;
  /** A figure that keeps its own height, centered with its caption under it. */
  fixed?: (paging: Paging) => ReactNode;
  /** The figure may continue onto further slides. */
  continues: boolean;
}

type Slide = { kind: "title"; id: "title"; title: string } | FigureSlide;

/** One step of the deck: a slide, or one page of a slide that continues. */
interface Frame {
  slide: Slide;
  page: number;
  pages: number;
}

type Stage = (part: SlidePart, node: ReactNode, grow: boolean) => ReactNode;

/** The beats a figure slide stages, in order: its answer, its figure, its
 *  caption. A continuation stages only its figure; the rest is already read. */
const partsOf = (s: FigureSlide): SlidePart[] => [...(s.result === undefined ? [] : (["result"] as const)), "figure", "caption"];

const frameTitle = (f: Frame): string => (f.pages > 1 ? `${f.slide.title} (${f.page + 1} of ${f.pages})` : f.slide.title);

const TitleSlide = ({ manifest, timeline }: { manifest: Manifest; timeline: SlideTimeline }) => {
  const audit = manifestAudit(manifest);
  return (
    <section aria-labelledby="slide-title" className="flex h-full min-h-0 flex-col items-start justify-center gap-6">
      <GrmMark size={112} motion="present" />
      <Reveal at={1} timeline={timeline}>
        <h2 id="slide-title" className="text-display font-semibold uppercase tracking-widest text-ink">
          Genomic Resistance Mapping
        </h2>
      </Reveal>
      <Reveal at={2} timeline={timeline}>
        <p className="type-headline">What GRM found</p>
        <p className="mt-2 flex flex-wrap items-center gap-x-2 gap-y-1 text-note text-muted-text">
          <span>Report built {dateText(manifest.generated_at)}</span>
          {audit.total === 0 ? null : (
            <>
              <span aria-hidden="true">·</span>
              {audit.invalid > 0 ? <strong className="font-semibold text-coral-text">Provisional metrics</strong> : <span>Internal cross-validation only</span>}
            </>
          )}
        </p>
      </Reveal>
    </section>
  );
};

/** Measures its box and hands the size to the figure, so a slide figure
 *  fills what the slide leaves and nothing scrolls. */
const FitBox = ({ figure, live }: { figure: (fit: Fit) => ReactNode; live: boolean }) => {
  const ref = useRef<HTMLDivElement | null>(null);
  const [fit, setFit] = useState<Fit | null>(null);
  useLayoutEffect(() => {
    const el = ref.current;
    if (el === null) return undefined;
    const read = () => setFit((f) => (f !== null && f.width === el.clientWidth && f.height === el.clientHeight ? f : { width: el.clientWidth, height: el.clientHeight }));
    read();
    const observer = new ResizeObserver(read);
    observer.observe(el);
    return () => observer.disconnect();
  }, []);
  return (
    <div ref={ref} data-figure-fit={live ? "" : undefined} className="relative min-h-0 min-w-0 flex-1">
      {fit !== null ? figure(fit) : null}
    </div>
  );
};

/** Off screen: whether a fixed figure and its caption fit the slide's box
 *  as one page; two pages when they do not and the figure marks itself
 *  data-continues. An unlaid box reports one. */
const FixedMeasure = ({ figure, caption, onPages }: { figure: ReactNode; caption: string; onPages(pages: number): void }) => {
  const box = useRef<HTMLDivElement | null>(null);
  const content = useRef<HTMLDivElement | null>(null);
  useLayoutEffect(() => {
    const b = box.current;
    const c = content.current;
    if (b === null || c === null) return undefined;
    const read = () => {
      const style = getComputedStyle(b);
      const room = b.clientHeight - (parseFloat(style.paddingTop) || 0) - (parseFloat(style.paddingBottom) || 0);
      const overflows = b.clientHeight > 0 && c.offsetHeight > room + 1;
      onPages(overflows && c.querySelector("[data-continues]") !== null ? 2 : 1);
    };
    read();
    const observer = new ResizeObserver(read);
    observer.observe(b);
    observer.observe(c);
    return () => observer.disconnect();
  }, [onPages]);
  return (
    <div ref={box} data-fixed-box="" className={FIXED_BOX}>
      <div ref={content} className="min-w-0">
        <div className="min-w-0">{figure}</div>
        <p className="prose-measure mt-6 text-note text-muted-text">{caption}</p>
      </div>
    </div>
  );
};

/** One figure slide's layout: title, answer, figure, caption. The live
 *  slide stages each beat; the measuring copy reports its page count. */
const SlideBody = ({ slide, paging, stage, measure }: { slide: FigureSlide; paging: Paging; stage: Stage; measure?: (pages: number) => void }) => {
  const live = measure === undefined;
  const Tag = live ? "section" : "div";
  const fill = slide.fill;
  return (
    <Tag aria-labelledby={live ? `${slide.id}-title` : undefined} className="flex h-full min-h-0 flex-col">
      <h2 id={live ? `${slide.id}-title` : undefined} className="text-body font-semibold text-muted-text">
        {slide.title}
      </h2>
      {/* A phone's slide sets the answer a size down, so the figure keeps
          the height its labels need; figure labels never shrink. */}
      {slide.result === undefined
        ? null
        : stage("result", <p className="mt-2 max-w-[56rem] text-[1.5rem]/[1.15] font-semibold tracking-[-0.014em] text-ink sm:text-stat">{slide.result}</p>, false)}
      {fill !== undefined ? (
        <>
          {stage(
            "figure",
            <div data-figure={live ? "" : undefined} className="mt-4 flex h-full min-h-0 flex-1 flex-col">
              <FitBox live={live} figure={(fit) => fill({ ...fit, page: paging.page, pages: paging.pages, onPages: measure })} />
            </div>,
            true,
          )}
          {stage(
            "caption",
            <p data-caption={live ? "" : undefined} className="prose-measure mt-3 text-note text-muted-text">
              {slide.caption}
            </p>,
            false,
          )}
        </>
      ) : measure !== undefined ? (
        <FixedMeasure figure={slide.fixed?.(paging)} caption={slide.caption} onPages={measure} />
      ) : (
        <div data-fixed-box="" className={FIXED_BOX}>
          {stage(
            "figure",
            <div data-figure className="min-w-0">
              {slide.fixed?.(paging)}
            </div>,
            false,
          )}
          {stage(
            "caption",
            <p data-caption className="prose-measure mt-6 text-note text-muted-text">
              {slide.caption}
            </p>,
            false,
          )}
        </div>
      )}
    </Tag>
  );
};

/** The measuring copy's beats: each in the same box the live slide's reveal
 *  would give it, shown at once. */
const plain: Stage = (_part, node, grow) => <div className={grow ? "flex min-h-0 flex-1 flex-col" : undefined}>{node}</div>;

const SlideView = ({
  frame,
  manifest,
  direction,
  timelineRef,
  onStage,
}: {
  frame: Frame;
  manifest: Manifest;
  direction: 1 | -1;
  timelineRef: MutableRefObject<SlideTimeline | null>;
  onStage(stage: number, steps: number): void;
}) => {
  const { slide } = frame;
  const parts = slide.kind === "figure" ? (frame.page > 0 ? (["figure"] as SlidePart[]) : partsOf(slide)) : [];
  // The title slide waits for the mark to assemble before the words land.
  const timeline = useSlideTimeline(slide.kind === "title" ? 2 : parts.length, slide.kind === "title" ? [900, 900] : undefined);
  useEffect(() => {
    timelineRef.current = timeline;
    onStage(timeline.stage, timeline.steps);
  });
  const stage: Stage = (part, node, grow) => {
    const at = parts.indexOf(part);
    const className = grow ? "flex min-h-0 flex-1 flex-col" : undefined;
    return at < 0 ? (
      <div className={className}>{node}</div>
    ) : (
      <Reveal at={at + 1} timeline={timeline} className={className}>
        {node}
      </Reveal>
    );
  };
  return (
    <div className={`h-full min-h-0 ${direction === 1 ? "slide-step-forward" : "slide-step-back"}`}>
      {slide.kind === "title" ? <TitleSlide manifest={manifest} timeline={timeline} /> : <SlideBody slide={slide} paging={{ page: frame.page, pages: frame.pages }} stage={stage} />}
    </div>
  );
};

export const PresentMode = ({ manifest, onClose }: { manifest: Manifest; onClose(): void }) => {
  const slides: Slide[] = useMemo(() => {
    const pick = ruleDataset(manifest);
    const lineage = pick === null ? null : ruleLineageScore(pick.row);
    return [
      { kind: "title", id: "title", title: "Genomic Resistance Mapping" },
      ...primerSteps(manifest).map(
        (s): Slide => ({ kind: "figure", id: s.id, title: s.title, result: s.lead, caption: s.caption, fixed: () => s.figure, continues: false }),
      ),
      ...FINDINGS.map((build): Slide => {
        const f = build(manifest, { page: false });
        const figure = f.figure;
        return {
          kind: "figure",
          id: f.id,
          title: f.title,
          result: f.result,
          caption: f.caption,
          ...(typeof figure === "function" ? { fill: figure } : { fixed: () => figure }),
          continues: f.continues === true,
        };
      }),
      {
        kind: "figure",
        id: "rule",
        title: "One interpretable rule",
        result:
          pick === null
            ? "No set covering machine rule in this report."
            : `The set covering machine reaches ${pct(pick.ba)} balanced accuracy for ${pick.row.antibiotic} in ${speciesShort(pick.row.species)} on random splits${lineage === null ? "" : `, ${pct(lineage)} with lineages held out`}.`,
        caption:
          "Carriage of the first rule's k-mers, by lab call. A rule marks k-mers that separate the calls in these isolates; they can travel with the cause on a plasmid or in a lineage rather than be the cause.",
        fixed: ({ page, pages }) => <RuleFigure manifest={manifest} page={page} pages={pages} />,
        continues: true,
      },
    ];
  }, [manifest]);

  // Pages per continuing slide, as the off-screen copies measure them.
  const [measured, setMeasured] = useState<Record<string, number>>({});
  const continuing = slides.filter((s): s is FigureSlide => s.kind === "figure" && s.continues);
  const continuingKey = continuing.map((s) => s.id).join(" ");
  const reporters = useMemo(
    () =>
      Object.fromEntries(
        continuingKey
          .split(" ")
          .filter((id) => id !== "")
          .map((id) => [id, (n: number) => setMeasured((m) => (m[id] === n ? m : { ...m, [id]: n }))]),
      ) as Record<string, (n: number) => void>,
    [continuingKey],
  );
  const frames: Frame[] = useMemo(
    () =>
      slides.flatMap((slide) => {
        const pages = slide.kind === "figure" && slide.continues ? Math.max(1, measured[slide.id] ?? 1) : 1;
        return Array.from({ length: pages }, (_, page) => ({ slide, page, pages }));
      }),
    [slides, measured],
  );

  // The place in the deck is a slide and a page, so a measurement that adds
  // or drops pages elsewhere never moves the reader.
  const [at, setAt] = useState<{ id: string; page: number }>({ id: "title", page: 0 });
  const [direction, setDirection] = useState<1 | -1>(1);
  const [beat, setBeat] = useState({ stage: 0, steps: 0 });
  const timelineRef = useRef<SlideTimeline | null>(null);
  const onStage = useCallback((stage: number, steps: number) => setBeat((b) => (b.stage === stage && b.steps === steps ? b : { stage, steps })), []);
  const found = frames.findIndex((f) => f.slide.id === at.id && f.page === Math.min(at.page, f.pages - 1));
  const index = Math.max(0, found);
  const frame = frames[index]!;
  const last = index === frames.length - 1;

  const jump = useCallback(
    (next: number) => {
      const clamped = Math.max(0, Math.min(next, frames.length - 1));
      const to = frames[clamped]!;
      setDirection(clamped >= index ? 1 : -1);
      setAt({ id: to.slide.id, page: to.page });
    },
    [frames, index],
  );

  const advance = useCallback(() => {
    if (timelineRef.current?.skipToEnd() === true) return;
    if (last) onClose();
    else jump(index + 1);
  }, [index, jump, last, onClose]);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (owns(e.target, KEY_OWNERS) || !(e.target instanceof Element) || e.target.closest("dialog")?.classList.contains("presentation-dialog") !== true) return;
      const act: Record<string, () => void> = {
        ArrowRight: advance,
        PageDown: advance,
        ArrowLeft: () => jump(index - 1),
        PageUp: () => jump(index - 1),
        Home: () => jump(0),
        End: () => jump(frames.length - 1),
      };
      const fn = act[e.key];
      if (fn !== undefined) {
        e.preventDefault();
        fn();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [advance, index, jump, frames.length]);

  return (
    <ModalDialog label="GRM findings presentation" onClose={onClose} className="presentation-dialog">
      <div className="flex h-full min-w-0 flex-col bg-paper text-ink">
        <header className={`${COLUMN} flex shrink-0 items-center justify-between gap-4 py-3`}>
          <span className="flex min-w-0 items-center gap-2.5">
            <GrmMark size={20} />
            <span data-deck-title="" className="min-w-0 break-words text-tiny font-semibold uppercase leading-snug tracking-widest text-blue">
              What GRM found · Genomic Resistance Mapping
            </span>
          </span>
          <button type="button" onClick={onClose} className={ACTION_BUTTON}>
            Close
          </button>
        </header>
        <main className={`${COLUMN} relative min-h-0 flex-1 py-4`}>
          <SlideView key={`${frame.slide.id}-${frame.page}`} frame={frame} manifest={manifest} direction={direction} timelineRef={timelineRef} onStage={onStage} />
          <div aria-hidden="true" inert className="pointer-events-none invisible absolute inset-0 overflow-hidden">
            {continuing.map((s) => (
              <div key={s.id} className="absolute inset-0 px-4 py-4 sm:px-8">
                <SlideBody slide={s} paging={{ page: 0, pages: undefined }} stage={plain} measure={reporters[s.id]} />
              </div>
            ))}
          </div>
        </main>
        <footer className="shrink-0 border-t border-hair">
          {/* One row from sm up; on a phone the outline takes the first row
              and Previous and Next the second, under the thumb. */}
          <div className={`${COLUMN} grid grid-cols-[auto_minmax(0,1fr)_auto] items-center gap-x-3 gap-y-2 py-3`}>
            <button type="button" onClick={() => jump(index - 1)} disabled={index === 0} className={`${ACTION_BUTTON} col-start-1 row-start-2 justify-self-start sm:row-start-1`}>
              Previous
            </button>
            <span className="col-span-3 col-start-1 row-start-1 flex min-w-0 items-center justify-center gap-3 sm:col-span-1 sm:col-start-2">
              <span role="status" aria-live="polite" className="sr-only">
                Slide {index + 1} of {frames.length}: {frameTitle(frame)}
              </span>
              {/* The control gives way first: its title wraps, the ticks
                  keep their width. */}
              <span className="min-w-0">
                <Dropdown
                  ariaLabel="Go to slide"
                  options={frames.map((f, i) => ({ id: String(i), label: `${i + 1} · ${frameTitle(f)}` }))}
                  mode={{ kind: "one", value: String(index), onSelect: (v) => jump(Number(v)) }}
                  wrap
                />
              </span>
              {beat.steps > 1 ? (
                <span aria-hidden="true" data-beats="" className="flex shrink-0 items-center gap-1">
                  {Array.from({ length: beat.steps }, (_, s) => (
                    <span key={s} className={`block size-1.5 rounded-full transition-colors duration-fast ease-out ${s < beat.stage ? "bg-blue" : "bg-ink/[0.14]"}`} />
                  ))}
                </span>
              ) : null}
            </span>
            <button type="button" onClick={advance} className={`${ACTION_BUTTON} col-start-3 row-start-2 justify-self-end sm:row-start-1`}>
              {last ? "Finish" : "Next"}
            </button>
          </div>
        </footer>
      </div>
    </ModalDialog>
  );
};
