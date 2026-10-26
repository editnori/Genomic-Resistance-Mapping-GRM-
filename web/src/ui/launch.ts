// Where the reader launched an overlay from: the centre of the row or mark
// last clicked or entered. The dataset inspector reads it once as it opens,
// so it grows from that row; a stale point is ignored.

interface Launch {
  x: number;
  y: number;
  at: number;
}

let last: Launch | null = null;

const FRESH_MS = 1000;

export const recordLaunch = (el: Element | null): void => {
  if (el === null) return;
  const r = el.getBoundingClientRect();
  last = { x: r.left + r.width / 2, y: r.top + r.height / 2, at: performance.now() };
};

export const takeLaunch = (): { x: number; y: number } | null => {
  const l = last;
  last = null;
  return l !== null && performance.now() - l.at < FRESH_MS ? { x: l.x, y: l.y } : null;
};
