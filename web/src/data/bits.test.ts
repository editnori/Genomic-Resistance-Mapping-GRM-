import { describe, expect, it } from "vitest";
import { decodeBits, encodeBits } from "./bits";

const ones = (bits: boolean[]) => bits.flatMap((b, i) => (b ? [i] : []));

describe("feature bits", () => {
  it("decodes the pipeline's own base64 test vector", () => {
    // src/analysis.rs: base64(&[0x0000_0000_0061_6263]) == "Y2JhAAAAAAA=".
    // Byte 0 is 0x63, byte 1 0x62, byte 2 0x61, little endian.
    expect(ones(decodeBits("Y2JhAAAAAAA=", 64))).toEqual([0, 1, 5, 6, 9, 13, 14, 16, 21, 22]);
  });

  it("reads bits past the first word in the right order", () => {
    const present = Array.from({ length: 130 }, (_, i) => i === 0 || i === 64 || i === 129);
    expect(ones(decodeBits(encodeBits(present), 130))).toEqual([0, 64, 129]);
  });

  it("ignores padding past the sample count", () => {
    expect(decodeBits(encodeBits([true, false, true]), 3)).toEqual([true, false, true]);
  });
});
