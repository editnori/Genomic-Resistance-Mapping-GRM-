import type { Config } from "tailwindcss";

/**
 * Theme colors flow through CSS variables set by `.theme-light` /
 * `.theme-dark` on <html> (src/styles/index.css).
 *
 * Channel-triplet vars (`--paper: 247 247 248`) are wrapped in
 * `rgb(var(..) / <alpha-value>)` so opacity modifiers work (`bg-blue/10`).
 * Full-value vars (`hair`, `bezel`, the fill ramp, the chart mark) carry
 * their own alpha and do not support the `/alpha` modifier.
 */
const channel = (name: string): string => `rgb(var(${name}) / <alpha-value>)`;

const INTER = ['"Inter Variable"', '"Inter"', "system-ui", "-apple-system", '"Segoe UI"', "sans-serif"];

// The lawful ranks, shared by their own names and the stock names below.
const DISPLAY: [string, { lineHeight: string; letterSpacing: string }] = ["2rem", { lineHeight: "1.12", letterSpacing: "-0.014em" }];
const STAT: [string, { lineHeight: string; letterSpacing: string }] = ["1.875rem", { lineHeight: "1.12", letterSpacing: "-0.014em" }];
const HEADLINE: [string, { lineHeight: string; letterSpacing: string }] = ["1.125rem", { lineHeight: "1.4", letterSpacing: "-0.014em" }];
const BODY: [string, { lineHeight: string }] = ["0.9375rem", { lineHeight: "1.5" }];
const NOTE: [string, { lineHeight: string }] = ["0.84375rem", { lineHeight: "1.5" }];
const CONTROL: [string, { lineHeight: string }] = ["0.8125rem", { lineHeight: "1.45" }];
const TINY: [string, { lineHeight: string }] = ["0.71875rem", { lineHeight: "1.35" }];

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
        bezel: "var(--bezel)",
        // One dataset or one learner score on a chart: the blue series.
        mark: "var(--chart-mark)",
      },
      fontFamily: { sans: INTER },
      // Fixed product ranks. Viewport width never changes the hierarchy.
      //   display  32px    view titles
      //   stat     30px    primary figures
      //   headline 18px    panel answers, one sentence
      //   body     15px    reading copy
      //   note     13.5px  supporting sentences
      //   control  13px    controls and table cells
      //   tiny     11.5px  short labels and state only
      // Stock size names resolve onto the nearest lawful rank.
      fontSize: {
        display: DISPLAY,
        stat: STAT,
        headline: HEADLINE,
        body: BODY,
        note: NOTE,
        control: CONTROL,
        tiny: TINY,
        xs: TINY,
        sm: NOTE,
        base: BODY,
        lg: HEADLINE,
        xl: HEADLINE,
        "2xl": STAT,
        "3xl": STAT,
        "4xl": DISPLAY,
        "5xl": DISPLAY,
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
