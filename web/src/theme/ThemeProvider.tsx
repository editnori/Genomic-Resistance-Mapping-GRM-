import { createContext, useCallback, useContext, useEffect, useMemo, useState, type ReactNode } from "react";
import { Moon, Sun } from "lucide-react";
import { FOCUS_RING } from "../ui/primitives";
import { THEME_STORAGE_KEY } from "./storage";

export type Theme = "light" | "dark";

interface ThemeContextValue {
  theme: Theme;
  toggleTheme(): void;
}

const ThemeContext = createContext<ThemeContextValue | null>(null);

function readStoredTheme(): Theme {
  // Light is the default; dark is an explicit reader choice.
  try {
    return localStorage.getItem(THEME_STORAGE_KEY) === "dark" ? "dark" : "light";
  } catch {
    return "light";
  }
}

export function ThemeProvider({ children }: { children: ReactNode }) {
  const [theme, setTheme] = useState<Theme>(readStoredTheme);

  useEffect(() => {
    const root = document.documentElement;
    root.classList.remove("theme-light", "theme-dark");
    root.classList.add(theme === "dark" ? "theme-dark" : "theme-light");
    root.style.colorScheme = theme;
    try {
      localStorage.setItem(THEME_STORAGE_KEY, theme);
    } catch {
      // Persistence is best effort; the in-memory theme still applies.
    }
  }, [theme]);

  const toggleTheme = useCallback(() => setTheme((t) => (t === "light" ? "dark" : "light")), []);
  const value = useMemo(() => ({ theme, toggleTheme }), [theme, toggleTheme]);
  return <ThemeContext.Provider value={value}>{children}</ThemeContext.Provider>;
}

export function useTheme(): ThemeContextValue {
  const ctx = useContext(ThemeContext);
  if (ctx === null) throw new Error("useTheme must be used within <ThemeProvider>");
  return ctx;
}

/** Squared and outlined like Present. The icon names the theme the button
 *  turns to; on a switch the outgoing icon turns away as the other turns in. */
export function ThemeToggle() {
  const { theme, toggleTheme } = useTheme();
  const dark = theme === "dark";
  const icon = (shown: boolean) => `theme-icon absolute h-4 w-4 ${shown ? "rotate-0 scale-100 opacity-100" : "-rotate-90 scale-50 opacity-0"}`;
  return (
    <button
      type="button"
      onClick={toggleTheme}
      aria-label={dark ? "Switch to light theme" : "Switch to dark theme"}
      className={`ui-control-target relative inline-flex h-8 w-8 items-center justify-center rounded-md border border-hair bg-paper-pure text-ink transition-[color,border-color,transform] duration-fast ease-out hover:border-blue hover:text-blue active:scale-[0.98] ${FOCUS_RING}`}
    >
      <Sun aria-hidden="true" className={icon(dark)} strokeWidth={1.75} />
      <Moon aria-hidden="true" className={icon(!dark)} strokeWidth={1.75} />
    </button>
  );
}
