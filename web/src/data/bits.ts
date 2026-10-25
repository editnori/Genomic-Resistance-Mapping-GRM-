// Feature presence arrives as base64 of little-endian u64 words where bit i
// is sample i. Little-endian words laid end to end are the same bits as a
// byte string where bit i sits in byte i >> 3 at position i & 7, so the
// decoded bytes are indexed directly.

export const decodeBase64 = (b64: string): Uint8Array => {
  const text = atob(b64);
  const out = new Uint8Array(text.length);
  for (let i = 0; i < text.length; i += 1) out[i] = text.charCodeAt(i);
  return out;
};

/** Presence over `n` samples as one boolean per sample. */
export const decodeBits = (b64: string, n: number): boolean[] => {
  const bytes = decodeBase64(b64);
  const out = new Array<boolean>(n);
  for (let i = 0; i < n; i += 1) out[i] = (((bytes[i >> 3] ?? 0) >> (i & 7)) & 1) === 1;
  return out;
};

/** Inverse of decodeBits, for fixtures and tests: whole u64 words, padded. */
export const encodeBits = (present: readonly boolean[]): string => {
  const words = Math.ceil(present.length / 64);
  const bytes = new Uint8Array(Math.max(words, 1) * 8);
  present.forEach((p, i) => {
    if (p) bytes[i >> 3]! |= 1 << (i & 7);
  });
  let text = "";
  bytes.forEach((b) => {
    text += String.fromCharCode(b);
  });
  return btoa(text);
};
