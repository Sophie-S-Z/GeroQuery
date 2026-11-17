import { Suspense, lazy, useEffect, useRef, useState } from "react";
import { flushSync } from "react-dom";
import GeneView from "./components/GeneView.jsx";
import Overview from "./components/Overview.jsx";
import { json } from "./lib/db.js";
import { fmtInt } from "./lib/format.js";
import { useUrlState } from "./lib/state.js";
import { useTheme } from "./lib/theme.js";

// The landing view and the gene view ship in the entry bundle: one is where a
// reader arrives, the other is the primary action and the target of the `/` and
// ⌘K shortcuts, which must find its search field the moment they fire. Every
// other view is fetched the first time it is opened.
const CrossSpecies = lazy(() => import("./components/CrossSpecies.jsx"));
const GeneList = lazy(() => import("./components/GeneList.jsx"));
const Landscape = lazy(() => import("./components/Landscape.jsx"));
const Mortality = lazy(() => import("./components/Mortality.jsx"));
const Panel = lazy(() => import("./components/Panel.jsx"));
const Perturbation = lazy(() => import("./components/Perturbation.jsx"));
const Volcano = lazy(() => import("./components/Volcano.jsx"));

/**
 * Apply a state change inside a view transition, so the outgoing view leaves
 * while the incoming one arrives and the tab marker slides between tabs.
 * Without browser support, or for a reader who asked for reduced motion, the
 * change simply happens.
 */
function withTransition(change) {
  const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  if (!document.startViewTransition || reduced) {
    change();
    return;
  }
  document.startViewTransition(() => flushSync(change));
}

const VIEWS = [
  ["overview", "Overview"],
  ["gene", "Gene"],
  ["list", "Gene list"],
  ["landscape", "Landscape"],
  ["certainty", "Certainty"],
  ["crossspecies", "Cross-species"],
  ["panel", "Panel"],
  ["mortality", "Mortality"],
  ["perturbation", "Perturbation"],
];

function Masthead({ meta, theme, onToggleTheme, onHome }) {
  return (
    <header className="masthead">
      <button type="button" className="wordmark" onClick={onHome}>
        GeroQuery
      </button>
      {meta ? (
        <div className="masthead-stats">
          <span>
            <b>{fmtInt(meta.n_contrasts_rows)}</b> effect sizes
          </span>
          <span>
            <b>{fmtInt(meta.n_genes)}</b> genes
          </span>
          <span>
            <b>{meta.n_studies}</b> contrasts / <b>{meta.n_series}</b> series
          </span>
          <span>
            <b>{meta.n_pinned_artifacts}</b> pinned artifacts
          </span>
        </div>
      ) : null}
      <span className="masthead-spacer" />
      {meta ? (
        <span
          className="mode-flag"
          data-mode={meta.signature_mode}
          title={
            meta.signature_mode === "full"
              ? "Built from the full signature table."
              : "Built from the committed curated slice, not the full table. Counts are smaller than the reported corpus."
          }
        >
          {meta.signature_mode === "full" ? "Full corpus" : "Curated slice"}
        </span>
      ) : null}
      <button type="button" className="theme-toggle" onClick={onToggleTheme}>
        {theme === "dark" ? "Light" : "Dark"}
      </button>
    </header>
  );
}

/** Shown for the moment a view's code is in flight: the shape, not a spinner. */
function ViewSkeleton() {
  return (
    <div className="panel" aria-busy="true" aria-label="Loading view">
      <div className="skeleton" style={{ width: "45%", height: "1.6rem" }} />
      <div className="skeleton" style={{ width: "70%" }} />
      <div className="skeleton" style={{ width: "60%" }} />
    </div>
  );
}

export default function App() {
  const [state, update] = useUrlState();
  const [meta, setMeta] = useState(null);
  const [, toggleTheme, resolvedTheme] = useTheme();

  useEffect(() => {
    json("meta.json").then(setMeta).catch(() => setMeta(null));
  }, []);

  // "/" focuses gene search from any view — the convention researchers already
  // know from GitHub and documentation sites — and is ignored while typing in a
  // field. ⌘K / Ctrl+K does the same from anywhere, as in Linear and Notion.
  const focusSearch = useRef(false);
  useEffect(() => {
    const onKey = (event) => {
      const commandK =
        event.key.toLowerCase() === "k" && (event.metaKey || event.ctrlKey) && !event.altKey;
      const slash = event.key === "/" && !event.metaKey && !event.ctrlKey && !event.altKey;
      if (!commandK && !slash) return;
      const { target } = event;
      const typing =
        target.isContentEditable || /^(INPUT|TEXTAREA|SELECT)$/.test(target.tagName);
      if (slash && typing) return;
      event.preventDefault();
      const field = document.querySelector("[data-search]");
      if (field) field.focus();
      else {
        focusSearch.current = true;
        update({ view: "gene" });
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [update]);

  useEffect(() => {
    if (!focusSearch.current) return;
    focusSearch.current = false;
    document.querySelector("[data-search]")?.focus();
  }, [state.view]);

  // Charts resolve CSS custom properties at draw time, so a theme change has to
  // force a redraw. Keying the view subtree on the theme is the cheapest
  // correct way to do that.
  const themeKey = resolvedTheme;

  return (
    <div className="app">
      <Masthead
        meta={meta}
        theme={resolvedTheme}
        onToggleTheme={() => withTransition(toggleTheme)}
        onHome={() => update({ view: "overview" })}
      />

      <nav className="tabs" role="tablist" aria-label="Views">
        {VIEWS.map(([key, label]) => (
          <button
            key={key}
            type="button"
            role="tab"
            className="tab"
            aria-selected={state.view === key}
            onClick={() => withTransition(() => update({ view: key }))}
          >
            {label}
            {state.view === key ? <span className="tab-indicator" aria-hidden="true" /> : null}
          </button>
        ))}
      </nav>

      <main key={themeKey}>
        <Suspense fallback={<ViewSkeleton />}>
          {state.view === "overview" ? (
            <Overview
              onNavigate={(view, patch) => update({ view, ...patch })}
              onPick={(gene) => update({ view: "gene", gene })}
            />
          ) : null}

          {state.view === "gene" ? (
            <GeneView
              symbol={state.gene}
              species={state.species}
              onPick={(gene) => update({ gene })}
              onSpecies={(species) => update({ species })}
            />
          ) : null}

          {state.view === "list" ? (
            <GeneList
              genes={state.genes}
              species={state.species}
              onUpdate={update}
              onPick={(gene) => update({ view: "gene", gene })}
            />
          ) : null}

          {state.view === "landscape" ? (
            <Landscape
              species={state.species}
              direction={state.dir}
              onPick={(gene) => update({ view: "gene", gene })}
              onSpecies={(species) => update({ species })}
              onDirection={(dir) => update({ dir })}
            />
          ) : null}

          {state.view === "perturbation" ? (
            <Perturbation onPick={(gene) => update({ view: "gene", gene })} />
          ) : null}

          {state.view === "certainty" ? (
            <Volcano
              species={state.species}
              onPick={(gene) => update({ view: "gene", gene })}
              onSpecies={(species) => update({ species })}
            />
          ) : null}

          {state.view === "crossspecies" ? (
            <CrossSpecies onPick={(gene) => update({ view: "gene", gene })} />
          ) : null}

          {state.view === "panel" ? <Panel /> : null}
          {state.view === "mortality" ? <Mortality /> : null}
        </Suspense>
      </main>

      <footer>
        Nothing here is asserted that was not measured. Every pooled effect is a random-effects
        meta-analysis of per-study standardized differences, computed in Python and read by this
        page — the browser renders numbers, it never derives them.
        {meta ? (
          <>
            {" "}
            Manifest {meta.manifest_version}, checksums verified {meta.checksums_verified_on}. Data
            version <span className="mono">{meta.data_version}</span>. Built {meta.built_on}.
          </>
        ) : null}
      </footer>
    </div>
  );
}
