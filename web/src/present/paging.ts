// When a slide takes a second page for its caption. The deck's off-screen
// copy of each figure slide lays it out as one page and measures four
// heights; this decides from them alone, so the rule is the same for every
// slide and a longer sentence on any slide cannot draw its figure over its
// caption. A figure squeezed below the height its own labels need draws at
// that least height and overflows its box, so the height it draws is what
// it needs. When the four do not fit, page 1 keeps the title, the answer
// and the figure, and page 2 shows the title and the caption: the figure
// stays beside the sentence it supports, and gets every pixel the caption
// gave up. A figure that marks itself data-continues moves part of itself
// to page 2 instead and keeps the caption on both pages.

/** One figure slide laid out as one page, in CSS pixels. */
export interface SlideFit {
  /** The slide's height. */
  room: number;
  /** From the slide's top to the figure's top: the title, the answer and
   *  the gap under them. */
  head: number;
  /** The figure as drawn in the box one page leaves it: its least height
   *  when that box is smaller. */
  figure: number;
  /** The caption with the gap above it and any space the slide keeps under
   *  it. */
  caption: number;
}

/** How a slide is paged: how many pages, and whether the caption is on
 *  each page or only on the last, alone under the title. */
export interface SlidePlan {
  pages: number;
  caption: "each" | "last";
}

/** A layout off by less than this is rounding, not overflow. */
export const SLACK = 1;

export const ONE_PAGE: SlidePlan = { pages: 1, caption: "each" };

/** Whether the title, the answer, the figure and the caption fit one page.
 *  An unlaid slide (no height yet) fits. */
export const fitsOnePage = ({ room, head, figure, caption }: SlideFit): boolean => room <= 0 || head + figure + caption <= room + SLACK;

/** The plan for a slide that does not page its own rows. `continues`: its
 *  figure marks itself data-continues and moves part of itself to page 2. */
export const slidePlan = (fit: SlideFit | null, continues: boolean): SlidePlan =>
  fit === null || fitsOnePage(fit) ? ONE_PAGE : continues ? { pages: 2, caption: "each" } : { pages: 2, caption: "last" };

export const samePlan = (a: SlidePlan | undefined, b: SlidePlan): boolean => a !== undefined && a.pages === b.pages && a.caption === b.caption;

/** What a slide's page shows: everything; the figure without the caption
 *  (a caption split's first page); or the title and the caption alone. */
export type PageShows = "all" | "figure" | "caption";

export const pageShows = (plan: SlidePlan, page: number): PageShows => (plan.caption === "each" ? "all" : page >= plan.pages - 1 ? "caption" : "figure");

const px = (v: string): number => parseFloat(v) || 0;

/** Reads a slide's one-page copy: `root` is the slide, `[data-measure-slot]`
 *  the figure's box and `[data-measure-caption]` the caption. The figure's
 *  extent is the lowest bottom of anything it draws, as a reader sees it
 *  (screen-reader text left out). Null until both parts are there. */
export const readSlideFit = (root: Element): SlideFit | null => {
  const slot = root.querySelector("[data-measure-slot]");
  const caption = root.querySelector("[data-measure-caption]");
  if (slot === null || caption === null) return null;
  const r = root.getBoundingClientRect();
  const s = slot.getBoundingClientRect();
  let bottom = s.top;
  for (const e of slot.querySelectorAll("*")) {
    if (e.classList.contains("sr-only") || e.closest('[role="tooltip"]') !== null) continue;
    const b = e.getBoundingClientRect();
    if (b.width > 0 || b.height > 0) bottom = Math.max(bottom, b.bottom);
  }
  // Space the slide keeps under the caption: the bottom padding of each box
  // between the caption and the slide.
  let under = 0;
  for (let e = caption.parentElement; e !== null && e !== root; e = e.parentElement) under += px(getComputedStyle(e).paddingBottom);
  return {
    room: r.height,
    head: s.top - r.top,
    figure: bottom - s.top,
    caption: caption.getBoundingClientRect().height + px(getComputedStyle(caption).marginTop) + under,
  };
};
