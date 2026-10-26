// Present: one deck on a full canvas. A title slide where the mark runs its
// full choreography, the five Overview findings, then the rule showcase.
// Each slide stages its beats (the answer, the figure, the caption) and the
// footer's ticks show them land; advancing mid-reveal completes the slide
// first. Left and right, Previous and Next, or a click on the slide move one
// step; the slide control opens the named outline for a direct jump. The
// header, slide, and footer share one column.

import { useCallback, useEffect, useRef, useState, type MutableRefObject, type ReactNode } from "react";
import type { Manifest } from "../data/schema";
import { dateText, pct, speciesShort } from "../data/format";
import { manifestAudit } from "../data/derive/evaluation";
import { GrmMark } from "../app/Mark";
import { ProvisionalChip } from "../ui/EvaluationNotice";
import { ModalDialog } from "../ui/ModalDialog";
import { Dropdown } from "../ui/Dropdown";
import { FigureBlock, PILL_ACTION, type FigureContent, type SlidePart } from "../ui/primitives";
import { FINDINGS } from "../views/findings";
import { Reveal, useSlideTimeline, type SlideTimeline } from "./motion";
import { RuleFigure, ruleDataset } from "./figures";

const COLUMN = "mx-auto w-full max-w-[72rem] px-4 sm:px-8";
const KEY_OWNERS = 'input, select, textarea, [role="combobox"], [role="listbox"], [role="option"], [contenteditable]:not([contenteditable="false"])';
/** A click on these does its own thing rather than advancing the slide. */
const CONTROLS = 'button, a, input, select, textarea, [role="combobox"], [role="listbox"], [role="option"], [role="grid"], [tabindex]';

const owns = (target: EventTarget | null, selector: string): boolean => target instanceof Element && target.closest(selector) !== null;

interface FigureSlide {
  kind: "figure";
  id: string;
  title: string;
  result?: string;
  figure: FigureContent;
  caption: string;
}

type Slide = { kind: "title"; id: "title"; title: string } | FigureSlide;

/** The beats a figure slide stages, in order: its answer, its figure, its
 *  caption. */
const partsOf = (s: FigureSlide): SlidePart[] => [...(s.result === undefined ? [] : (["result"] as const)), "figure", "caption"];

const TitleSlide = ({ manifest, timeline }: { manifest: Manifest; timeline: SlideTimeline }) => (
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
        <ProvisionalChip audit={manifestAudit(manifest)} />
        <span>Report built {dateText(manifest.generated_at)}</span>
      </p>
    </Reveal>
  </section>
);

const SlideView = ({
  slide,
  manifest,
  direction,
  timelineRef,
  onStage,
}: {
  slide: Slide;
  manifest: Manifest;
  direction: 1 | -1;
  timelineRef: MutableRefObject<SlideTimeline | null>;
  onStage(stage: number, steps: number): void;
}) => {
  const parts = slide.kind === "figure" ? partsOf(slide) : [];
  // The title slide waits for the mark to assemble before the words land.
  const timeline = useSlideTimeline(slide.kind === "title" ? 2 : parts.length, slide.kind === "title" ? [900, 900] : undefined);
  useEffect(() => {
    timelineRef.current = timeline;
    onStage(timeline.stage, timeline.steps);
  });
  const stage = (part: SlidePart, node: ReactNode, grow: boolean) => (
    <Reveal at={parts.indexOf(part) + 1} timeline={timeline} className={grow ? "flex min-h-0 flex-1 flex-col" : undefined}>
      {node}
    </Reveal>
  );
  return (
    <div className={`h-full min-h-0 ${direction === 1 ? "slide-step-forward" : "slide-step-back"}`}>
      {slide.kind === "title" ? (
        <TitleSlide manifest={manifest} timeline={timeline} />
      ) : (
        <FigureBlock id={slide.id} title={slide.title} result={slide.result} figure={slide.figure} caption={slide.caption} variant="slide" stage={stage} />
      )}
    </div>
  );
};

export const PresentMode = ({ manifest, onClose }: { manifest: Manifest; onClose(): void }) => {
  const pick = ruleDataset(manifest);
  const slides: Slide[] = [
    { kind: "title", id: "title", title: "Genomic Resistance Mapping" },
    ...FINDINGS.map((finding): Slide => ({ kind: "figure", ...finding(manifest, { page: false }) })),
    {
      kind: "figure",
      id: "rule",
      title: "One interpretable rule",
      result: pick === null ? "No set covering machine rule in this report." : `The set covering machine reaches ${pct(pick.ba)} balanced accuracy for ${pick.row.antibiotic} in ${speciesShort(pick.row.species)}.`,
      figure: <RuleFigure manifest={manifest} />,
      caption: "Carriage of the first rule's k-mers, by lab call. Balanced accuracy is the whole model's, out of fold on random splits.",
    },
  ];
  const [index, setIndex] = useState(0);
  const [direction, setDirection] = useState<1 | -1>(1);
  const [beat, setBeat] = useState({ stage: 0, steps: 0 });
  const timelineRef = useRef<SlideTimeline | null>(null);
  const onStage = useCallback((stage: number, steps: number) => setBeat((b) => (b.stage === stage && b.steps === steps ? b : { stage, steps })), []);
  const slide = slides[index]!;
  const last = index === slides.length - 1;

  const jump = useCallback(
    (next: number) => {
      const clamped = Math.max(0, Math.min(next, slides.length - 1));
      setDirection(clamped >= index ? 1 : -1);
      setIndex(clamped);
    },
    [index, slides.length],
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
        End: () => jump(slides.length - 1),
      };
      const fn = act[e.key];
      if (fn !== undefined) {
        e.preventDefault();
        fn();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [advance, index, jump, slides.length]);

  return (
    <ModalDialog label="GRM findings presentation" onClose={onClose} className="presentation-dialog">
      <div className="flex h-full min-w-0 flex-col bg-paper text-ink">
        <header className={`${COLUMN} flex shrink-0 items-center justify-between gap-4 py-3`}>
          <span className="flex min-w-0 items-center gap-2.5">
            <GrmMark size={20} />
            <span className="truncate text-tiny font-semibold uppercase tracking-widest text-blue">What GRM found · Genomic Resistance Mapping</span>
          </span>
          <button type="button" onClick={onClose} className={PILL_ACTION}>
            Close
          </button>
        </header>
        <main
          className={`${COLUMN} min-h-0 flex-1 py-4`}
          onClick={(e) => {
            if (!owns(e.target, CONTROLS)) advance();
          }}
        >
          <SlideView key={slide.id} slide={slide} manifest={manifest} direction={direction} timelineRef={timelineRef} onStage={onStage} />
        </main>
        <footer className="shrink-0 border-t border-hair">
          <div className={`${COLUMN} flex items-center justify-between gap-3 py-3`}>
            <button type="button" onClick={() => jump(index - 1)} disabled={index === 0} className={PILL_ACTION}>
              Previous
            </button>
            <span className="flex min-w-0 items-center gap-3">
              <span role="status" aria-live="polite" className="sr-only">
                Slide {index + 1} of {slides.length}: {slide.title}
              </span>
              {/* The control gives way first: its title truncates, the
                  ticks keep their width. */}
              <span className="min-w-0">
                <Dropdown
                  ariaLabel="Go to slide"
                  options={slides.map((s, i) => ({ id: String(i), label: `${i + 1} · ${s.title}` }))}
                  mode={{ kind: "one", value: String(index), onSelect: (v) => jump(Number(v)) }}
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
            <button type="button" onClick={advance} className={PILL_ACTION}>
              {last ? "Finish" : "Next"}
            </button>
          </div>
        </footer>
      </div>
    </ModalDialog>
  );
};
