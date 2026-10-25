import type { Config } from "tailwindcss";

/**
 * Theme colors flow through CSS variables set by `.theme-light` /
 * `.theme-dark` on <html> (src/styles/index.css).
 *
 * Channel-triplet vars (`--paper: 247 247 248`) are wrapped in
 * `rgb(var(..) / <alpha-value>)` so opacity modifiers work (`bg-blue/10`).
 * Full-value vars (`hair`, the fill ramp, the status pair) carry their own
 * alpha and do not support the `/alpha` modifier.
 */
const channel = (name: string): string => `rgb(var(${name}) / <alpha-value>)`;

const INTER = ['"Inter Variable"', '"Inter"', "system-ui", "-apple-system", '"Segoe UI"', "sans-serif"];

export default {
  content: ["./index.html", "./src/**/*.{ts,tsx}"],
  theme: {
    extend: {
      colors: {
        paper: channel("--paper"),
        "paper-pure": channel("--paper-pure"),
        "sidebar-surface": channel("--sidebar-surface"),
        "surface-subtle": "var(--fill-1)",
        "surface-control": "var(--fill-2)",
        "surface-strong": "var(--fill-3)",
        ink: channel("--ink"),
        blue: {
          DEFAULT: channel("--blue"),
          body: channel("--blue-body"),
          deep: channel("--blue-deep"),
        },
        coral: channel("--coral"),
        live: channel("--live"),
        // Text-safe accent channels: the same hue families as the coral/live
        // marks, tuned per theme so small text stays readable on its paper.
        "live-text": channel("--live-text"),
        "coral-text": channel("--coral-text"),
        "muted-text": channel("--muted-text"),
        hair: "var(--hair)",
        // A verdict against a stated bar: solid green meets, quiet magenta
        // misses. Never a phenotype and never an interaction color.
        meet: { DEFAULT: "var(--status-meet-bg)", fg: "var(--status-meet-fg)" },
        miss: { DEFAULT: "var(--status-miss-bg)", fg: "var(--status-miss-fg)" },
        mark: "var(--chart-mark)",
      },
      fontFamily: { sans: INTER },
      // Fixed product ranks, size/line height. Viewport width never changes
      // the hierarchy.
      //   display  32/36  page titles
      //   stat     30/34  a Present result
      //   headline 18/25  a figure's result sentence
      //   body     15/22  figure titles and reading copy
      //   note     13/20  captions
      //   control  13/18  tables and controls
      //   tiny     12/16  axes and secondary values
      // Stock size names resolve onto the nearest lawful rank.
      fontSize: {
        display: ["2rem", { lineHeight: "2.25rem", letterSpacing: "-0.014em" }],
        stat: ["1.875rem", { lineHeight: "2.125rem", letterSpacing: "-0.014em" }],
        headline: ["1.125rem", { lineHeight: "1.5625rem", letterSpacing: "-0.011em" }],
        body: ["0.9375rem", { lineHeight: "1.375rem" }],
        note: ["0.8125rem", { lineHeight: "1.25rem" }],
        control: ["0.8125rem", { lineHeight: "1.125rem" }],
        tiny: ["0.75rem", { lineHeight: "1rem" }],
        xs: ["0.75rem", { lineHeight: "1rem" }],
        sm: ["0.8125rem", { lineHeight: "1.25rem" }],
        base: ["0.9375rem", { lineHeight: "1.375rem" }],
        lg: ["1.125rem", { lineHeight: "1.5625rem", letterSpacing: "-0.011em" }],
        xl: ["1.125rem", { lineHeight: "1.5625rem", letterSpacing: "-0.011em" }],
        "2xl": ["1.875rem", { lineHeight: "2.125rem", letterSpacing: "-0.014em" }],
        "3xl": ["1.875rem", { lineHeight: "2.125rem", letterSpacing: "-0.014em" }],
        "4xl": ["2rem", { lineHeight: "2.25rem", letterSpacing: "-0.014em" }],
        "5xl": ["2rem", { lineHeight: "2.25rem", letterSpacing: "-0.014em" }],
      },
      // Weight ceiling 600: the heavy utilities keep their names and clamp.
      fontWeight: { bold: "600", extrabold: "600", black: "600" },
      letterSpacing: { tight: "var(--track-title)", tighter: "var(--track-title)" },
      // Motion rides the CSS clocks; a bare `transition-*` uses the fast clock.
      transitionDuration: {
        DEFAULT: "var(--dur-fast)",
        fast: "var(--dur-fast)",
        med: "var(--dur-med)",
        slow: "var(--dur-slow)",
      },
      transitionTimingFunction: { DEFAULT: "var(--ease-out)", out: "var(--ease-out)" },
    },
  },
  plugins: [],
} satisfies Config;
