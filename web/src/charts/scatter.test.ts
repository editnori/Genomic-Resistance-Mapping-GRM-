import { describe, expect, it } from "vitest";
import { coreBounds, flatMapNote } from "./ScatterMap";

describe("coreBounds", () => {
  it("frames the main cluster and leaves a distant one outside", () => {
    // 90 genomes spread over [0, 0.09], 10 genomes in a cluster at 1.
    const xs = [...Array.from({ length: 90 }, (_, i) => (i % 10) * 0.01), ...Array.from({ length: 10 }, () => 1)];
    const ys = [...Array.from({ length: 90 }, (_, i) => Math.floor(i / 10) * 0.01), ...Array.from({ length: 10 }, () => 1)];
    const b = coreBounds(xs, ys)!;
    const inside = xs.filter((x, i) => x >= b[0] && x <= b[1] && ys[i]! >= b[2] && ys[i]! <= b[3]).length;
    expect(inside).toBe(90);
  });

  it("gives no frame when every genome sits on one point", () => {
    expect(coreBounds([0, 0, 0], [1, 1, 1])).toBeNull();
  });
});

describe("flatMapNote", () => {
  it("names the narrow vertical axis of a whole-map frame and gives the eigenvalue share", () => {
    // Spans 2 and 0.03; sums of squares 2 and 0.0006.
    const xs = [-1, 0, 1];
    const ys = [0.01, -0.02, 0.01];
    expect(flatMapNote(xs, ys, xs, ys)).toBe(
      "The coordinates themselves form a line: drawn to one scale, this frame's vertical span is 1.5% of its horizontal span. Over all genomes, the second principal coordinate's eigenvalue, the sum of its squared coordinates, is 0.03% of the first's.",
    );
  });

  it("names the horizontal axis when a zoomed frame is the narrow one, and keeps the eigenvalue over all genomes", () => {
    // A distant genome at -3.0003 carries the first axis; the frame holds the
    // other three, which share almost one horizontal position. Sums of squares
    // over all genomes: about 12 and 0.18.
    const xs = [-3.0003, 1, 1.0003, 1];
    const ys = [0, -0.3, 0, 0.3];
    expect(flatMapNote(xs.slice(1), ys.slice(1), xs, ys)).toBe(
      "The coordinates themselves form a line: drawn to one scale, this frame's horizontal span is 0.05% of its vertical span. Over all genomes, the second principal coordinate's eigenvalue, the sum of its squared coordinates, is 1.5% of the first's.",
    );
  });

  it("says nothing when both axes of the frame are visible, at the 5% floor or above", () => {
    expect(flatMapNote([-1, 0, 1], [-0.05, 0, 0.05], [-1, 0, 1], [-0.05, 0, 0.05])).toBeNull();
    expect(flatMapNote([-1, 0, 1], [-1, 1, 0], [-1, 0, 1], [-1, 1, 0])).toBeNull();
  });

  it("says nothing for fewer than three genomes in the frame or a frame on one point", () => {
    expect(flatMapNote([0, 1], [0, 0], [0, 1, 2], [0, 0, 1])).toBeNull();
    expect(flatMapNote([1, 1, 1], [2, 2, 2], [0, 1, 1, 1], [0, 2, 2, 2])).toBeNull();
  });
});
