import { THEME_STORAGE_KEY } from "./theme/storage";

try {
  if (localStorage.getItem(THEME_STORAGE_KEY) === "dark") {
    document.documentElement.classList.replace("theme-light", "theme-dark");
    document.documentElement.style.colorScheme = "dark";
  }
} catch {
  // Storage can be unavailable in hardened browser contexts; light stays.
}
