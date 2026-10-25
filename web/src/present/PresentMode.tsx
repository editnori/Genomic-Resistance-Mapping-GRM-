// One six-slide deck: the five Overview findings, then the rule showcase.
import { useEffect, useState } from "react";
import type { Manifest } from "../data/schema";
import { int, pct, speciesShort } from "../data/format";
import { manifestAudit } from "../data/derive/evaluation";
import { Brand } from "../app/Brand";
import { ReportDetails } from "../app/ReportDetails";
import { ModalDialog } from "../ui/ModalDialog";
import { FigureBlock, FOCUS_RING } from "../ui/primitives";
import { FINDINGS } from "../views/findings";
import { RuleFigure, ruleDataset } from "./figures";

const CHROME = `inline-flex min-h-11 items-center justify-center rounded-md px-3 text-control font-medium text-ink hover:bg-surface-control disabled:cursor-not-allowed disabled:opacity-45 ${FOCUS_RING}`;
const KEY_OWNERS = 'input, select, textarea, [role="combobox"], [role="listbox"], [role="option"], [contenteditable]:not([contenteditable="false"])';

export const PresentMode = ({ manifest, onClose }: { manifest: Manifest; onClose(): void }) => {
  const [index, setIndex] = useState(0);
  const [direction, setDirection] = useState<1 | -1>(1);
  const pick = ruleDataset(manifest);
  const slides = [
    ...FINDINGS.map((finding) => finding(manifest, { page: false })),
    {
      id: "rule",
      title: "One interpretable rule",
      result: pick === null ? "No set covering machine rule in this report." : `The set covering machine reaches ${pct(pick.ba)} balanced accuracy for ${pick.row.antibiotic} in ${speciesShort(pick.row.species)}.`,
      figure: <RuleFigure manifest={manifest} />,
      caption: "The first rule's k-mer carriage, by lab call. Balanced accuracy is the whole model's, out of fold on random splits. Carriage is association, not proof of mechanism.",
    },
  ];
  const slide = slides[index]!;
  const audit = manifestAudit(manifest);
  const jump = (next: number) => {
    const clamped = Math.max(0, Math.min(next, slides.length - 1));
    setDirection(clamped >= index ? 1 : -1);
    setIndex(clamped);
  };
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (!(e.target instanceof Element) || e.target.closest(KEY_OWNERS) !== null || !e.target.closest("dialog")?.classList.contains("presentation-dialog")) return;
      const next = e.key === "ArrowRight" || e.key === "PageDown" ? index + 1
        : e.key === "ArrowLeft" || e.key === "PageUp" ? index - 1
          : e.key === "Home" ? 0 : e.key === "End" ? slides.length - 1 : null;
      if (next === null) return;
      e.preventDefault();
      setDirection(next >= index ? 1 : -1);
      setIndex(Math.max(0, Math.min(next, slides.length - 1)));
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [index, slides.length]);
  return (
    <ModalDialog label="GRM findings presentation" onClose={onClose} className="presentation-dialog">
      <div className="flex h-full min-w-0 flex-col bg-paper text-ink">
        <header className="flex shrink-0 items-center justify-between gap-4 px-4 py-2 sm:px-8">
          <Brand />
          <div className="flex items-center gap-1"><ReportDetails manifest={manifest} /><button type="button" onClick={onClose} className={CHROME}>Close</button></div>
        </header>
        <main className="mx-auto min-h-0 w-full max-w-[72rem] flex-1 px-4 py-4 sm:px-8">
          <div key={slide.id} className={`h-full min-h-0 ${direction === 1 ? "slide-step-forward" : "slide-step-back"}`}>
            <FigureBlock {...slide} variant="slide" />
          </div>
        </main>
        <footer className="shrink-0 border-t border-hair px-4 pb-2 sm:px-8">
          <div className="flex items-center justify-between gap-2">
            <button type="button" onClick={() => jump(index - 1)} disabled={index === 0} className={CHROME}>Previous</button>
            <label className="flex min-w-0 items-center gap-2 text-control">
              <span className="sr-only">Go to slide</span>
              <select value={index} onChange={(e) => jump(Number(e.target.value))} className={`min-h-11 max-w-[11rem] truncate rounded-md bg-paper px-1 text-ink ${FOCUS_RING}`}>
                {slides.map((s, i) => <option key={s.id} value={i}>{i + 1} / {slides.length} · {s.title}</option>)}
              </select>
            </label>
            <button type="button" onClick={() => index === slides.length - 1 ? onClose() : jump(index + 1)} className={CHROME}>{index === slides.length - 1 ? "Finish" : "Next"}</button>
          </div>
          <p className="text-center text-tiny text-muted-text">Report {manifest.generated_at.slice(0, 10)}{manifest.evaluation_revision === undefined ? "" : ` · revision ${manifest.evaluation_revision}`} · {audit.total === 0 ? "No evaluations" : audit.invalid > 0 ? `${int(audit.invalid)} of ${int(audit.total)} provisional` : "Internal CV only"}</p>
          <span role="status" aria-live="polite" className="sr-only">Slide {index + 1} of {slides.length}: {slide.title}</span>
        </footer>
      </div>
    </ModalDialog>
  );
};
