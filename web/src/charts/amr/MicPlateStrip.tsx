// Bench-first broth microdilution plate strip in two modes:
// - mode="cutoff" (primer): shows doubling wells and the vertical S|R cutoff
//   line alone, with no essential-agreement band.
// - mode="agreement" (Overview MIC and Dataset Inspector): marks one
//   illustrative measured well and highlights its two immediate ±1 doubling
//   neighbours, apart from any S|R cutoff line.

import { mgL } from "../../data/format";

export interface MicPlateStripProps {
  /** "cutoff" draws the S|R call boundary alone (no ±1 tolerance band);
   *  "agreement" marks one measured well and its two ±1 doubling neighbours. */
  mode?: "cutoff" | "agreement";
  /** Resistant-above breakpoint as a log2 dilution step (matching
   *  `breakpoint.r_above` in the dataset report, where 0 = 1 mg/L and -1 = 0.5 mg/L).
   *  Printed unrounded via `mgL(rAbove)`. Omit when folds use different thresholds. */
  rAbove?: number | null;
  /** Label printed above the cutoff divider when `rAbove` is given, e.g. "Cohort breakpoint". */
  cutoffLabel?: string;
  /** Compact mode hides the outer caption when embedded inside a tight slot. */
  compact?: boolean;
  /** Optional short label above the plate. */
  label?: string;
}

export const MicPlateStrip = ({ mode = "agreement", rAbove, cutoffLabel = "Cohort breakpoint", compact = false, label }: MicPlateStripProps) => {
  const isReal = rAbove !== undefined && rAbove !== null && Number.isFinite(rAbove);
  const cutLabel = isReal ? mgL(rAbove) : null;

  if (mode === "cutoff") {
    const baseStep = isReal ? Math.floor(rAbove) : -1;
    const wellSteps = [baseStep - 2, baseStep - 1, baseStep, baseStep + 1, baseStep + 2, baseStep + 3];
    const wellLabels = wellSteps.map((step) => mgL(step));
    return (
      <figure data-amr-primitive="mic-plate" data-plate-mode="cutoff" className="m-0 w-full max-w-[36rem]">
        {label !== undefined ? <p className="mb-1.5 text-tiny font-medium text-muted-text">{label}</p> : null}
        <svg
          viewBox="0 0 640 138"
          role="img"
          aria-label={
            isReal
              ? `Drawing of a six-well broth microdilution plate (${wellLabels[0]} to ${wellLabels[5]} mg/L) with the susceptible and resistant cutoff above ${cutLabel} mg/L.`
              : "Drawing, not data: six-well broth microdilution plate (0.125 to 4 mg/L) showing the vertical cutoff that separates susceptible calls on the left from resistant calls on the right."
          }
          className="block h-auto w-full"
        >
          <rect x="4" y="24" width="632" height="66" rx="10" fill="var(--fill-1)" stroke="var(--chart-axis)" strokeWidth="1.25" />
          <text x="320" y="15" textAnchor="middle" fontSize="12" fontWeight="600" fill="rgb(var(--ink))">
            {isReal ? `${cutoffLabel} > ${cutLabel} mg/L` : "S | R cutoff"}
          </text>
          <line x1="320" y1="20" x2="320" y2="94" stroke="rgb(var(--ink))" strokeWidth="2" />
          {wellLabels.map((wellText, i) => {
            const cx = 64 + i * 102.4;
            const isS = i === 2;
            const isR = i === 3;
            const lit = isS || isR;
            const fill = isS ? "var(--series-s)" : isR ? "var(--series-r)" : "rgb(var(--paper-pure))";
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
                  {isS ? "S" : isR ? "R" : wellText}
                </text>
                {lit ? (
                  <text x={cx} y="108" textAnchor="middle" fontSize="11.5" fontWeight="500" fill="rgb(var(--muted-text))">
                    {wellText} mg/L
                  </text>
                ) : null}
              </g>
            );
          })}
          <text x="166" y="131" textAnchor="middle" fontSize="12" fontWeight="500" fill="rgb(var(--muted-text))">
            Susceptible (S)
          </text>
          <text x="474" y="131" textAnchor="middle" fontSize="12" fontWeight="500" fill="rgb(var(--muted-text))">
            Resistant (R)
          </text>
        </svg>
        {compact ? null : (
          <figcaption className="mt-1.5 text-tiny text-muted-text">
            Drawing, not data: labs test doubling drug levels (mg/L) and call the isolate susceptible (S) or resistant (R) against the cutoff.
          </figcaption>
        )}
      </figure>
    );
  }

  // Agreement mode: mark one illustrative measured well and its two immediate
  // ±1 doubling neighbours, completely apart from any S|R cutoff line.
  const baseStep = isReal ? Math.floor(rAbove) : -1;
  // When a real cohort breakpoint is shown, place it between wells 1 and 2
  // (x = 217.6) and place the illustrative measured well at index 4 with its
  // ±1 neighbours at 3 and 5, well clear of the cutoff line. Without a cutoff,
  // center the measured well at index 3 with ±1 neighbours at 2 and 4.
  const wellSteps = isReal
    ? [baseStep - 1, baseStep, baseStep + 1, baseStep + 2, baseStep + 3, baseStep + 4]
    : [-3, -2, -1, 0, 1, 2];
  const wellLabels = wellSteps.map((step) => mgL(step));
  const measuredIdx = isReal ? 4 : 3;
  const bandStartX = 64 + (measuredIdx - 1) * 102.4 - 36;
  const bandWidth = 2 * 102.4 + 72;
  const measuredCx = 64 + measuredIdx * 102.4;

  return (
    <figure data-amr-primitive="mic-plate" data-plate-mode="agreement" className="m-0 w-full max-w-[36rem]">
      {label !== undefined ? <p className="mb-1.5 text-tiny font-medium text-muted-text">{label}</p> : null}
      <svg
        viewBox="0 0 640 140"
        role="img"
        aria-label={
          isReal
            ? `Drawing of a six-well broth microdilution plate (${wellLabels[0]} to ${wellLabels[5]} mg/L) with ${cutoffLabel.toLowerCase()} above ${cutLabel} mg/L; well ${wellLabels[measuredIdx]} mg/L is marked as an illustrative measured MIC with its ±1 doubling neighbours (${wellLabels[measuredIdx - 1]} and ${wellLabels[measuredIdx + 1]} mg/L) highlighted apart from the cutoff.`
            : `Drawing, not data: six-well broth microdilution plate (${wellLabels[0]} to ${wellLabels[5]} mg/L) with ${wellLabels[measuredIdx]} mg/L marked as an illustrative measured MIC and its two immediate neighbours (${wellLabels[measuredIdx - 1]} and ${wellLabels[measuredIdx + 1]} mg/L) highlighted as ±1 doubling well essential agreement.`
        }
        className="block h-auto w-full"
      >
        <rect x="4" y="24" width="632" height="66" rx="10" fill="var(--fill-1)" stroke="var(--chart-axis)" strokeWidth="1.25" />

        {/* Optional cohort breakpoint divider, well apart from the ±1 band */}
        {isReal ? (
          <>
            <text x="166" y="15" textAnchor="middle" fontSize="11.5" fontWeight="600" fill="rgb(var(--ink))">
              {`${cutoffLabel} > ${cutLabel} mg/L`}
            </text>
            <line x1="217.6" y1="20" x2="217.6" y2="94" stroke="rgb(var(--ink))" strokeWidth="2" />
          </>
        ) : null}

        {/* ±1 doubling well essential-agreement band around one measured MIC */}
        <rect
          x={bandStartX}
          y="30"
          width={bandWidth}
          height="54"
          rx="8"
          fill="rgb(var(--blue) / 0.12)"
          stroke="rgb(var(--blue) / 0.65)"
          strokeWidth="1.5"
          strokeDasharray="5 3"
        />
        <text x={measuredCx} y="15" textAnchor="middle" fontSize="11.5" fontWeight="600" fill="rgb(var(--blue))">
          Illustrative measured MIC
        </text>

        {/* Six doubling dilution wells */}
        {wellLabels.map((wellText, i) => {
          const cx = 64 + i * 102.4;
          const isMeasured = i === measuredIdx;
          const isNeighbour = Math.abs(i - measuredIdx) === 1;
          const fill = isMeasured ? "rgb(var(--blue))" : "rgb(var(--paper-pure))";
          const stroke = isMeasured ? "rgb(var(--blue))" : isNeighbour ? "rgb(var(--blue) / 0.7)" : "var(--chart-axis)";
          return (
            <g key={i}>
              <circle cx={cx} cy="57" r={isMeasured ? 20 : 18} fill={fill} stroke={stroke} strokeWidth={isMeasured ? 2 : 1.5} />
              <text
                x={cx}
                y="61"
                textAnchor="middle"
                fontSize="12"
                fontWeight={isMeasured ? "600" : "500"}
                fill={isMeasured ? "rgb(var(--paper-pure))" : "rgb(var(--ink))"}
              >
                {wellText}
              </text>
              {isMeasured ? (
                <text x={cx} y="109" textAnchor="middle" fontSize="11.5" fontWeight="500" fill="rgb(var(--muted-text))">
                  {wellText} mg/L
                </text>
              ) : null}
            </g>
          );
        })}

        {/* Bottom zone label */}
        {isReal ? (
          <text x="115" y="132" textAnchor="middle" fontSize="12" fontWeight="500" fill="rgb(var(--muted-text))">
            Susceptible
          </text>
        ) : null}
        <text x={measuredCx} y="132" textAnchor="middle" fontSize="12" fontWeight="600" fill="rgb(var(--blue))">
          ±1 doubling well (agrees)
        </text>
      </svg>
      {compact ? null : (
        <figcaption className="mt-1.5 text-tiny text-muted-text">
          {isReal
            ? `Drawing, not data: ${cutoffLabel.toLowerCase()} (> ${cutLabel} mg/L) shown on doubling steps; dashed box illustrates ±1 doubling well around one measured MIC, wherever it sits on the plate.`
            : "Drawing, not data: labs test doubling drug levels (mg/L); essential agreement counts a predicted MIC within one doubling step either side of an isolate's measured MIC, wherever it sits on the plate."}
        </figcaption>
      )}
    </figure>
  );
};
