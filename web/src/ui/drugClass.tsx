// Drug-class identity: the three largest antibiotic families each get one
// colour and one shape; every other class shares a neutral ring. The class
// name is always printed beside the mark, so the mark is a scanning aid and
// never the only way to read a class.
//
// Three is the most identity colours that pass the dataviz palette checks
// for every pair at once (a table row can sit beside any other), in both
// themes, on paper and paper-pure: worst all-pairs CVD ΔE 13.0 light / 13.2
// dark, normal-vision ΔE 19.6 / 19.3. A fourth passing set always needs
// green, which this UI already uses for "meets the bar". In light mode the
// yellow and magenta marks sit below 3:1 on paper (2.0 and 2.5:1), so the
// printed name is the required relief. Forced colors draw every mark in
// CanvasText, where the shape still tells them apart.

export type ClassFamily = "beta-lactam" | "aminoglycoside" | "quinolone" | "other";

export interface FamilyInfo {
  id: ClassFamily;
  /** The family's name in a key: "β-lactams". */
  label: string;
  /** The mark's shape: the channel that survives without colour. */
  shape: "circle" | "square" | "triangle" | "ring";
}

export const FAMILIES: readonly FamilyInfo[] = [
  { id: "beta-lactam", label: "β-lactams", shape: "circle" },
  { id: "aminoglycoside", label: "Aminoglycosides", shape: "square" },
  { id: "quinolone", label: "Quinolones", shape: "triangle" },
  { id: "other", label: "Other classes", shape: "ring" },
];

const BY_ID = new Map(FAMILIES.map((f) => [f.id, f]));

/** The family a curated class name belongs to. Penicillins, cephalosporins
 *  (any generation, with or without an inhibitor), carbapenems and
 *  monobactams are β-lactams; an empty, unknown or rarer class is other. */
export const classFamily = (drugClass: string | null | undefined): ClassFamily => {
  const c = (drugClass ?? "").toLowerCase();
  if (/penicillin|cephalosporin|carbapenem|monobactam|penem\b|beta-lactam|β-lactam/.test(c)) return "beta-lactam";
  if (/aminoglycoside/.test(c)) return "aminoglycoside";
  if (/quinolone/.test(c)) return "quinolone";
  return "other";
};

export const familyInfo = (id: ClassFamily): FamilyInfo => BY_ID.get(id)!;

/** Classes ordered for a list: by family in key order, then by the given
 *  count, most first, then by name. */
export const byFamily = (classes: readonly string[], count: (c: string) => number = () => 0): string[] => {
  const rank = (c: string) => FAMILIES.findIndex((f) => f.id === classFamily(c));
  return [...classes].sort((a, b) => rank(a) - rank(b) || count(b) - count(a) || a.localeCompare(b));
};

/** The family's mark: a 10px shape in the family's colour. Decorative; the
 *  class or family name beside it carries the meaning. */
export const ClassGlyph = ({ drugClass, family, size = 10, className = "" }: { drugClass?: string; family?: ClassFamily; size?: number; className?: string }) => {
  const f = family ?? classFamily(drugClass);
  const shape = familyInfo(f).shape;
  const fill = `var(--class-${f})`;
  return (
    <svg aria-hidden="true" width={size} height={size} viewBox="0 0 10 10" data-class-family={f} data-shape={shape} className={`class-glyph inline-block shrink-0 ${className}`}>
      {shape === "circle" ? <circle cx={5} cy={5} r={4.25} fill={fill} /> : null}
      {shape === "square" ? <rect x={1} y={1} width={8} height={8} rx={1} fill={fill} /> : null}
      {shape === "triangle" ? <path d="M5 0.75 9.4 9 0.6 9Z" fill={fill} /> : null}
      {shape === "ring" ? <circle cx={5} cy={5} r={3.5} fill="none" stroke={fill} strokeWidth={1.5} /> : null}
    </svg>
  );
};

/** A class name with its family mark before it, in the text's own ink.
 *  Inline, so a long name wraps and a line clamp still cuts it; the mark
 *  and the name's first word never break apart. */
export const ClassName = ({ drugClass, label, className = "" }: { drugClass: string; label?: string; className?: string }) => {
  const text = drugClass === "" ? "Class not recorded" : (label ?? drugClass);
  const cut = text.indexOf(" ");
  const first = cut === -1 ? text : text.slice(0, cut);
  return (
    <span className={className}>
      <span className="whitespace-nowrap">
        <ClassGlyph drugClass={drugClass} className="mr-1 align-[-0.04em]" />
        {first}
      </span>
      {cut === -1 ? null : text.slice(cut)}
    </span>
  );
};

/** What each mark means, printed once: the shape as it appears, then the
 *  family in words. */
export const ClassKey = ({ className = "" }: { className?: string }) => (
  <span data-class-key="" className={`inline-flex flex-wrap items-baseline gap-x-3 gap-y-1 ${className}`}>
    {FAMILIES.map((f) => (
      <span key={f.id} className="inline-flex items-baseline gap-1">
        <ClassGlyph family={f.id} className="self-center" />
        {f.label}
      </span>
    ))}
  </span>
);
