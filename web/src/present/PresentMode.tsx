// Present: one deck on a full canvas. A title slide where the mark runs its
// full choreography, the five primer steps (one real report read from zero),
// the nine Overview findings, then the rule showcase.
// Each slide stages its beats (the answer, the figure, the caption) and the
// footer's ticks show them land; advancing mid-reveal completes the slide
// first. Only the deck's own controls move it: Previous and Next, the arrow
// and page keys, and the outline for a direct jump; a click on a slide does
// nothing. A figure too dense for the canvas at its minimum row height
// continues onto further slides that repeat its title, answer, and caption,
// so nothing scrolls and no label shrinks. Every other figure slide whose
// title, answer, figure at its least height and caption do not fit one page
// takes a second page: page 1 keeps the title, the answer and the figure,
// and page 2 shows the title and the caption (paging.ts decides). An
// off-screen copy of each figure slide measures what it needs. The header,
// slide, and footer share one column. A fixed figure sits one set gap under
// the answer at every width rather than centred in what is left, so a short
// slide has no band of empty space above its figure. No title or outline
// label is cut off: on a phone they wrap, and the footer puts the outline on
// its own row above Previous and Next.

import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState, type MutableRefObject, type ReactNode } from "react";
import type { Manifest } from "../data/schema";
import { datasetTitle, dateText, pct } from "../data/format";
import { manifestAudit } from "../data/derive/evaluation";
import { GrmMark } from "../app/Mark";
import { MeasuringCopy } from "../charts/kit";
import { ModalDialog } from "../ui/ModalDialog";
import { Dropdown } from "../ui/Dropdown";
import { ACTION_BUTTON, type Fit, type SlidePart } from "../ui/primitives";
import { AUDIT_SCOPE } from "../ui/EvaluationNotice";
import { FINDINGS, type PagedFit } from "../views/findings";
import { Reveal, useSlideTimeline, type SlideTimeline } from "./motion";
import { RuleFigure, ruleDataset, ruleLineageScore } from "./figures";
import { primerSteps } from "../views/primer";
import { useMicBaselines } from "../data/load";
import { micConstants } from "../data/baselines";
import { ONE_PAGE, pageShows, readSlideFit, samePlan, slidePlan, type PageShows, type SlidePlan } from "./paging";

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
  /** A filling figure that pages its own rows (it reports its page count
   *  through onPages); a fixed figure says so with data-continues instead. */
  continues: boolean;
}

type Slide = { kind: "title"; id: "title"; title: string } | FigureSlide;

/** One step of the deck: a slide, or one page of a slide that continues. */
interface Frame {
  slide: Slide;
  page: number;
  pages: number;
  /** What this page shows of its slide. */
  shows: PageShows;
}

type Stage = (part: SlidePart, node: ReactNode, grow: boolean) => ReactNode;

/** The beats a figure slide's page stages, in order: its answer, its
 *  figure, its caption. A continuation stages only its figure; the rest is
 *  already read. A caption split's first page stages the answer and the
 *  figure, and its second page the caption alone. */
const partsOf = (s: FigureSlide, page: number, shows: PageShows): SlidePart[] => {
  if (shows === "caption") return ["caption"];
  if (page > 0) return ["figure"];
  return [...(s.result === undefined ? [] : (["result"] as const)), "figure", ...(shows === "all" ? (["caption"] as const) : [])];
};

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
        <p className="type-body mt-1">How well does a bacterium&apos;s genome predict whether a lab finds it resistant or susceptible to an antibiotic?</p>
        <p className="mt-2 flex flex-wrap items-center gap-x-2 gap-y-1 text-note text-muted-text">
          <span>Report built {dateText(manifest.generated_at)}</span>
          {audit.total === 0 ? null : (
            <>
              <span aria-hidden="true">·</span>
              {audit.invalid > 0 ? <strong className="font-semibold text-coral-text">Provisional metrics</strong> : <span>{AUDIT_SCOPE}</span>}
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
    <div ref={ref} data-figure-fit={live ? "" : undefined} data-measure-slot={live ? undefined : ""} className="relative min-h-0 min-w-0 flex-1">
      {fit !== null ? figure(fit) : null}
    </div>
  );
};

/** One figure slide's layout: title, answer, figure, caption. The live
 *  slide stages each beat and shows what its page holds; the measuring copy
 *  (`measuring`) lays out the whole slide as one page and marks the figure's
 *  box and the caption for readSlideFit. */
const SlideBody = ({
  slide,
  paging,
  stage,
  shows = "all",
  measuring = false,
  onPages,
}: {
  slide: FigureSlide;
  paging: Paging;
  stage: Stage;
  shows?: PageShows;
  measuring?: boolean;
  /** A measuring copy of a figure that pages its own rows: where it reports
   *  its page count. */
  onPages?(pages: number): void;
}) => {
  const live = !measuring;
  const Tag = live ? "section" : "div";
  const fill = slide.fill;
  const caption = (gap: string) =>
    stage(
      "caption",
      <p data-caption={live ? "" : undefined} data-measure-caption={live ? undefined : ""} className={`prose-measure ${gap} text-note text-muted-text`}>
        {slide.caption}
      </p>,
      false,
    );
  return (
    <Tag aria-labelledby={live ? `${slide.id}-title` : undefined} className="flex h-full min-h-0 flex-col">
      <h2 id={live ? `${slide.id}-title` : undefined} className="text-body font-semibold text-muted-text">
        {slide.title}
      </h2>
      {/* A phone's slide sets the answer a size down, so the figure keeps
          the height its labels need; figure labels never shrink. */}
      {slide.result === undefined || shows === "caption"
        ? null
        : stage("result", <p className="mt-2 max-w-[56rem] text-[1.5rem]/[1.15] font-semibold tracking-[-0.014em] text-ink sm:text-stat">{slide.result}</p>, false)}
      {shows === "caption" ? (
        caption("mt-3")
      ) : fill !== undefined ? (
        <>
          {stage(
            "figure",
            <div data-figure={live ? "" : undefined} className="mt-4 flex h-full min-h-0 flex-1 flex-col">
              <FitBox live={live} figure={(fit) => fill({ ...fit, page: paging.page, pages: paging.pages, onPages })} />
            </div>,
            true,
          )}
          {shows === "all" ? caption("mt-3") : null}
        </>
      ) : (
        <div data-fixed-box="" className={FIXED_BOX}>
          {stage(
            "figure",
            <div data-figure={live ? "" : undefined} data-measure-slot={live ? undefined : ""} className="min-w-0">
              {slide.fixed?.(paging)}
            </div>,
            false,
          )}
          {shows === "all" ? caption("mt-6") : null}
        </div>
      )}
    </Tag>
  );
};

/** The measuring copy's beats: each in the same box the live slide's reveal
 *  would give it, shown at once. */
const plain: Stage = (_part, node, grow) => <div className={grow ? "flex min-h-0 flex-1 flex-col" : undefined}>{node}</div>;

/** What the footer shows on a page: its outline label, its beat count, and
 *  whether it is the last page (Finish for Next). */
interface FooterFace {
  label: string;
  steps: number;
  last: boolean;
}

/** Off screen: one figure slide laid out as one page, above a copy of the
 *  footer row its first page shows, so the copy has the room that page has
 *  on screen whichever page is showing now. A figure that pages its own
 *  rows reports its page count itself; for every other slide the copy reads
 *  its four heights and slidePlan decides, again whenever the slide, its
 *  text or its figure changes size. */
const MeasuredSlide = ({ slide, footer, onPlan }: { slide: FigureSlide; footer: FooterFace; onPlan(plan: SlidePlan): void }) => {
  const box = useRef<HTMLDivElement | null>(null);
  const own = slide.fill !== undefined && slide.continues;
  const onPages = useCallback((pages: number) => onPlan({ pages, caption: "each" }), [onPlan]);
  useLayoutEffect(() => {
    const root = box.current?.firstElementChild;
    if (own || root === null || root === undefined) return undefined;
    const read = () => onPlan(slidePlan(readSlideFit(root), root.querySelector("[data-continues]") !== null));
    const resize = new ResizeObserver(read);
    const watch = () => {
      resize.disconnect();
      for (const e of [root, ...root.querySelectorAll("h2, p, [data-measure-slot], [data-measure-slot] > *")]) resize.observe(e);
    };
    // A figure that draws once its box is measured adds its parts later, and
    // one that redraws in a box of the same size changes their attributes.
    const added = new MutationObserver(() => {
      watch();
      read();
    });
    added.observe(root, { childList: true, subtree: true, attributes: true, characterData: true });
    watch();
    read();
    return () => {
      added.disconnect();
      resize.disconnect();
    };
  }, [own, onPlan]);
  return (
    <div className="absolute inset-0 flex flex-col">
      <div className={`${COLUMN} relative min-h-0 flex-1`}>
        <div ref={box} data-slide-measure={slide.id} className="absolute inset-0 px-4 py-4 sm:px-8">
          <MeasuringCopy.Provider value>
            <SlideBody slide={slide} paging={{ page: 0, pages: undefined }} stage={plain} measuring onPages={own ? onPages : undefined} />
          </MeasuringCopy.Provider>
        </div>
      </div>
      <div data-footer-measure="" className="shrink-0 border-t border-hair">
        <FooterRow options={[{ id: "0", label: footer.label }]} value={0} steps={footer.steps} stage={0} first={false} last={footer.last} />
      </div>
    </div>
  );
};

/** The footer's row: Previous, the outline with the beat ticks, Next. One
 *  row from sm up; on a phone the outline takes the first row and Previous
 *  and Next the second, under the thumb. The live footer and its off-screen
 *  copies share it, so a copy is as tall as the footer on its page. */
const FooterRow = ({
  options,
  value,
  steps,
  stage,
  first,
  last,
  status,
  onPrevious,
  onNext,
  onJump,
}: {
  options: { id: string; label: string }[];
  value: number;
  steps: number;
  stage: number;
  first: boolean;
  last: boolean;
  /** The live footer's words for a screen reader; a copy has none. */
  status?: string;
  onPrevious?(): void;
  onNext?(): void;
  onJump?(index: number): void;
}) => (
  <div className={`${COLUMN} grid grid-cols-[auto_minmax(0,1fr)_auto] items-center gap-x-3 gap-y-2 py-3`}>
    <button type="button" onClick={onPrevious} disabled={first} className={`${ACTION_BUTTON} col-start-1 row-start-2 justify-self-start sm:row-start-1`}>
      Previous
    </button>
    <span className="col-span-3 col-start-1 row-start-1 flex min-w-0 items-center justify-center gap-3 sm:col-span-1 sm:col-start-2">
      {status === undefined ? null : (
        <span role="status" aria-live="polite" className="sr-only">
          {status}
        </span>
      )}
      {/* The control gives way first: its title wraps, the ticks keep their
          width. */}
      <span className="min-w-0">
        <Dropdown
          ariaLabel="Go to slide"
          options={options}
          mode={{ kind: "one", value: String(value), onSelect: (v) => onJump?.(Number(v)) }}
          wrap
        />
      </span>
      {steps > 1 ? (
        <span aria-hidden="true" data-beats={status === undefined ? undefined : ""} className="flex shrink-0 items-center gap-1">
          {Array.from({ length: steps }, (_, s) => (
            <span key={s} className={`block size-1.5 rounded-full transition-colors duration-fast ease-out ${s < stage ? "bg-blue" : "bg-ink/[0.14]"}`} />
          ))}
        </span>
      ) : null}
    </span>
    <button type="button" onClick={onNext} className={`${ACTION_BUTTON} col-start-3 row-start-2 justify-self-end sm:row-start-1`}>
      {last ? "Finish" : "Next"}
    </button>
  </div>
);


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
  const parts = slide.kind === "figure" ? partsOf(slide, frame.page, frame.shows) : [];
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
      {slide.kind === "title" ? (
        <TitleSlide manifest={manifest} timeline={timeline} />
      ) : (
        // A caption split's figure is one page of its own.
        <SlideBody slide={slide} paging={frame.shows === "all" ? { page: frame.page, pages: frame.pages } : { page: 0, pages: 1 }} stage={stage} shows={frame.shows} />
      )}
    </div>
  );
};

export const PresentMode = ({ manifest, onClose }: { manifest: Manifest; onClose(): void }) => {
  const files = useMicBaselines(manifest);
  const constants = useMemo(() => (files.status === "ready" ? micConstants(manifest.datasets, files.data, manifest.evaluation_revision) : undefined), [files, manifest]);
  const slides: Slide[] = useMemo(() => {
    const pick = ruleDataset(manifest);
    const lineage = pick === null ? null : ruleLineageScore(pick.row);
    return [
      { kind: "title", id: "title", title: "Genomic Resistance Mapping" },
      ...primerSteps(manifest).map(
        (s): Slide => ({ kind: "figure", id: s.id, title: s.title, result: s.lead, caption: s.caption, fixed: () => s.figure, continues: false }),
      ),
      ...FINDINGS.map((build): Slide => {
        const f = build(manifest, { page: false, micConstants: constants });
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
        title: "One rule a person can read",
        result:
          pick === null
            ? "No set covering machine rule in this report."
            : `The set covering machine reaches ${pct(pick.ba)} balanced accuracy for ${datasetTitle(pick.row)} on random splits${lineage === null ? "" : `, ${pct(lineage)} with lineages held out`}.`,
        caption:
          "This learner writes short yes-or-no rules. The bars count isolates carrying the first rule's k-mers, by lab call. A k-mer can travel with the cause, on a plasmid or in a lineage, without being it.",
        fixed: ({ page, pages }) => <RuleFigure manifest={manifest} page={page} pages={pages} />,
        continues: true,
      },
    ];
  }, [manifest, constants]);

  // Each figure slide's plan, as its off-screen copy measures it.
  const [measured, setMeasured] = useState<Record<string, SlidePlan>>({});
  const figures = slides.filter((s): s is FigureSlide => s.kind === "figure");
  const figuresKey = figures.map((s) => s.id).join(" ");
  const reporters = useMemo(
    () =>
      Object.fromEntries(
        figuresKey
          .split(" ")
          .filter((id) => id !== "")
          .map((id) => [id, (plan: SlidePlan) => setMeasured((m) => (samePlan(m[id], plan) ? m : { ...m, [id]: plan }))]),
      ) as Record<string, (plan: SlidePlan) => void>,
    [figuresKey],
  );
  const frames: Frame[] = useMemo(
    () =>
      slides.flatMap((slide) => {
        const plan = slide.kind === "figure" ? (measured[slide.id] ?? ONE_PAGE) : ONE_PAGE;
        const pages = Math.max(1, plan.pages);
        return Array.from({ length: pages }, (_, page) => ({ slide, page, pages, shows: pageShows(plan, page) }));
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
  const outline = frames.map((f, i) => ({ id: String(i), label: `${i + 1} · ${frameTitle(f)}` }));
  // The footer a slide's first page shows. A figure that pages its own rows
  // is measured under its first page's label; any other slide under the
  // label it has as one page, the page its copy lays out.
  const footerOf = (s: FigureSlide): FooterFace => {
    const at = Math.max(0, frames.findIndex((f) => f.slide.id === s.id));
    const f = frames[at]!;
    return s.fill !== undefined && s.continues
      ? { label: outline[at]!.label, steps: partsOf(s, 0, f.shows).length, last: at === frames.length - 1 }
      : { label: `${at + 1} · ${s.title}`, steps: partsOf(s, 0, "all").length, last: s === slides[slides.length - 1] };
  };

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
        {/* The slide and the footer share this box, so the copies below can
            lay out over both. */}
        <div className="relative flex min-h-0 flex-1 flex-col">
          <main className={`${COLUMN} min-h-0 flex-1 py-4`}>
            <SlideView key={`${frame.slide.id}-${frame.page}`} frame={frame} manifest={manifest} direction={direction} timelineRef={timelineRef} onStage={onStage} />
            {/* Off screen, over the slide and the footer: each figure slide as
                one page, above the footer row its first page shows. */}
            <div aria-hidden="true" inert className="pointer-events-none invisible absolute inset-0 overflow-hidden">
              {figures.map((s) => (
                <MeasuredSlide key={s.id} slide={s} footer={footerOf(s)} onPlan={reporters[s.id]!} />
              ))}
            </div>
          </main>
          <footer className="shrink-0 border-t border-hair">
            <FooterRow
              options={outline}
              value={index}
              steps={beat.steps}
              stage={beat.stage}
              first={index === 0}
              last={last}
              status={`Slide ${index + 1} of ${frames.length}: ${frameTitle(frame)}`}
              onPrevious={() => jump(index - 1)}
              onNext={advance}
              onJump={jump}
            />
          </footer>
        </div>
      </div>
    </ModalDialog>
  );
};
