/**
 * Theme plumbing.
 *
 * Observable Plot draws into an SVG with resolved colour values, so it cannot
 * read CSS custom properties the way the rest of the page does. `css()` reads
 * the computed values off :root once per render, which keeps a single source of
 * truth in styles.css instead of a second palette hard-coded in JavaScript.
 */

import { useCallback, useEffect, useState } from "react";

const STORAGE_KEY = "geroquery-theme";

/** Resolved design tokens, for canvas/SVG contexts that cannot use var(). */
export function css() {
  const style = getComputedStyle(document.documentElement);
  const read = (name) => style.getPropertyValue(name).trim();
  return {
    ink: read("--ink"),
    inkSoft: read("--ink-soft"),
    inkFaint: read("--ink-faint"),
    paper: read("--paper"),
    paperRaised: read("--paper-raised"),
    rule: read("--rule"),
    up: read("--up"),
    down: read("--down"),
    // `null` is a reserved word, so the token is spelled with a trailing
    // underscore rather than quoted at every call site.
    null_: read("--null"),
    mono: read("--mono"),
  };
}

/**
 * Light/dark with an explicit override that beats the system preference in
 * both directions.
 *
 * @returns {[("light"|"dark"|"system"), () => void, ("light"|"dark")]}
 */
export function useTheme() {
  const [choice, setChoice] = useState(
    () => localStorage.getItem(STORAGE_KEY) ?? "system",
  );
  const [systemDark, setSystemDark] = useState(
    () => window.matchMedia("(prefers-color-scheme: dark)").matches,
  );

  useEffect(() => {
    const media = window.matchMedia("(prefers-color-scheme: dark)");
    const onChange = (event) => setSystemDark(event.matches);
    media.addEventListener("change", onChange);
    return () => media.removeEventListener("change", onChange);
  }, []);

  const resolved = choice === "system" ? (systemDark ? "dark" : "light") : choice;

  useEffect(() => {
    if (choice === "system") {
      document.documentElement.removeAttribute("data-theme");
      localStorage.removeItem(STORAGE_KEY);
    } else {
      document.documentElement.setAttribute("data-theme", choice);
      localStorage.setItem(STORAGE_KEY, choice);
    }
  }, [choice]);

  const toggle = useCallback(() => {
    setChoice(resolved === "dark" ? "light" : "dark");
  }, [resolved]);

  return [choice, toggle, resolved];
}
