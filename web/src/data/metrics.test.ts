import { describe, expect, it } from "vitest";
import { agreement, barVerdict, markerCalls, REFERENCE_BARS } from "./metrics";

describe("agreement", () => {
  it("counts VME over resistant isolates and ME over susceptible ones", () => {
    // 4 resistant (1 missed), 6 susceptible (2 called resistant).
    const truth = [true, true, true, true, false, false, false, false, false, false];
    const calls = [true, true, true, false, true, true, false, false, false, false];
    const a = agreement(truth, calls);
    expect([a.tp, a.fn, a.fp, a.tn]).toEqual([3, 1, 2, 4]);
    expect(a.vme).toBeCloseTo(1 / 4);
    expect(a.me).toBeCloseTo(2 / 6);
    expect(a.balanced_accuracy).toBeCloseTo((3 / 4 + 4 / 6) / 2);
  });
});

describe("markerCalls", () => {
  const a = [true, true, false, false];
  const b = [true, false, true, false];

  it("calls resistant when any marker is present", () => {
    expect(markerCalls([a, b], 4, "any")).toEqual([true, true, true, false]);
  });

  it("calls resistant only when every marker is present", () => {
    expect(markerCalls([a, b], 4, "all")).toEqual([true, false, false, false]);
  });

  it("calls nothing resistant with no markers chosen", () => {
    expect(markerCalls([], 3, "all")).toEqual([false, false, false]);
  });
});

describe("reference bars", () => {
  it("meets a bar at exactly the bar value and misses above it", () => {
    expect(barVerdict(0.015, REFERENCE_BARS.vme)).toBe("meets");
    expect(barVerdict(0.0151, REFERENCE_BARS.vme)).toBe("misses");
  });
});
