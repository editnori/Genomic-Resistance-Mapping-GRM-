import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import "@fontsource-variable/inter/opsz.css";
import "@fontsource-variable/inter/opsz-italic.css";
import "./styles/index.css";
import { ThemeProvider } from "./theme/ThemeProvider";
import { App } from "./app/App";
import { redirectMovedPath } from "./app/router";

const root = document.getElementById("root");
if (root === null) throw new Error("index.html is missing #root");

redirectMovedPath();

createRoot(root).render(
  <StrictMode>
    <ThemeProvider>
      <App />
    </ThemeProvider>
  </StrictMode>,
);
