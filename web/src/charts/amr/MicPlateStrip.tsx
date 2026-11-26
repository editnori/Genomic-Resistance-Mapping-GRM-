// Bench-first broth microdilution plate strip: shows doubling antibiotic
// dilutions across six wells, the vertical S|R cutoff line, and the ±1 well
// broth-dilution tolerance band. Parameterized with a dataset's actual
// breakpoint when available, or illustrative doubling steps otherwise.

import { mgL } from "../../data/format";

export interface MicPlateStripProps {
  /** Resistant-above breakpoint as a log2 dilution step (matching
   *  `breakpoint.r_above` in the dataset report, where 0 = 1 mg/L and -1 = 0.5 mg/L).
   *  When omitted, shows an illustrative 0.125–4 mg/L plate (cutoff > 0.5 mg/L). */
  rAbove?: number | null;
  /** Compact mode hides the outer caption when embedded inside a tight slot. */
  compact?: boolean;
  /** Optional short label above the plate. */
  label?: string;
}

/** Six doubling wells centered around the S | R cutoff in log2 steps:
 *  index 0,1 = clear susceptible; index 2 = S cutoff well (cutLog2);
 *  index 3 = R cutoff well (cutLog2 + 1); index 4,5 = clear resistant. */
export const MicPlateStrip = ({ rAbove, compact = false, label }: MicPlateStripProps) => {
  const isReal = rAbove !== undefined && rAbove !== null && Number.isFinite(rAbove);
  const cutLog2 = isReal ? Math.round(rAbove) : -1;
  const wellSteps = [cutLog2 - 2, cutLog2 - 1, cutLog2, cutLog2 + 1, cutLog2 + 2, cutLog2 + 3];
  const wellLabels = wellSteps.map((step) => mgL(step));
  const cutLabel = mgL(cutLog2);

  return (
    <figure data-amr-primitive="mic-plate" className="m-0 w-full max-w-[36rem]">
      {label !== undefined ? <p className="mb-1.5 text-tiny font-medium text-muted-text">{label}</p> : null}
      <svg
        viewBox="0 0 640 138"
        role="img"
        aria-label={
          isReal
            ? `Broth microdilution plate with six doubling steps from ${wellLabels[0]} to ${wellLabels[5]} mg/L; cutoff above ${cutLabel} mg/L separating susceptible and resistant, with the two cutoff wells highlighted as the ±1 dilution tolerance zone.`
            : "Drawing of a six-well broth microdilution plate (0.125 to 4 mg/L): two clear susceptible wells, two cutoff wells (±1 doubling dilution across the S and R boundary), and two clear resistant wells."
        }
        className="block h-auto w-full"
      >
        {/* Plate tray */}
        <rect x="4" y="24" width="632" height="66" rx="10" fill="var(--fill-1)" stroke="var(--chart-axis)" strokeWidth="1.25" />

        {/* ±1 dilution band around the S | R cutoff */}
        <rect
          x="216"
          y="30"
          width="208"
          height="54"
          rx="8"
          fill="rgb(var(--coral) / 0.14)"
          stroke="rgb(var(--coral-text) / 0.65)"
          strokeWidth="1.5"
          strokeDasharray="5 3"
        />

        {/* Cutoff divider */}
        <text x="320" y="15" textAnchor="middle" fontSize="12" fontWeight="600" fill="rgb(var(--ink))">
          {isReal ? `Cutoff > ${cutLabel} mg/L` : "Cutoff"}
        </text>
        <line x1="320" y1="20" x2="320" y2="94" stroke="rgb(var(--ink))" strokeWidth="2" />

        {/* Six doubling dilution wells */}
        {wellLabels.map((wellText, i) => {
          const cx = 64 + i * 102.4;
          const isCutoffS = i === 2;
          const isCutoffR = i === 3;
          const lit = isCutoffS || isCutoffR;
          const fill = isCutoffS ? "var(--series-s)" : isCutoffR ? "var(--series-r)" : "rgb(var(--paper-pure))";
          const stroke = lit ? fill : "var(--chart-axis)";
          return (
            <g key={i}>
              <circle cx={cx} cy="57" r={lit ? 20 : 18} fill={fill} stroke={stroke} strokeWidth={lit ? 2 : 1.5} />
              <text
                x={cx}
                y="61"
                textAnchor="middle"
                fontSize="12"
                fontWeight={lit ? "600" : "500"}
                fill={lit ? "rgb(var(--paper-pure))" : "rgb(var(--ink))"}
              >
                {isCutoffS ? "S" : isCutoffR ? "R" : wellText}
              </text>
              {lit ? (
                <text x={cx} y="108" textAnchor="middle" fontSize="11.5" fontWeight="500" fill="rgb(var(--muted-text))">
                  {wellText} mg/L
                </text>
              ) : null}
            </g>
          );
        })}

        {/* Bottom zone labels on dedicated row, well clear of cutoff numbers */}
        <text x="115" y="131" textAnchor="middle" fontSize="12" fontWeight="500" fill="rgb(var(--muted-text))">
          Susceptible
        </text>
        <text x="320" y="131" textAnchor="middle" fontSize="12" fontWeight="600" fill="rgb(var(--coral-text))">
          ±1 doubling well
        </text>
        <text x="525" y="131" textAnchor="middle" fontSize="12" fontWeight="500" fill="rgb(var(--muted-text))">
          Resistant
        </text>
      </svg>
      {compact ? null : (
        <figcaption className="mt-1.5 text-tiny text-muted-text">
          {isReal
            ? `Doubling steps in mg/L around the dataset breakpoint (${cutLabel} mg/L); essential agreement counts predictions within one step.`
            : "Drawing, not data: labs test doubling drug levels (mg/L); one doubling step either side of the measured well counts as essential agreement."}
        </figcaption>
      )}
    </figure>
  );
};
