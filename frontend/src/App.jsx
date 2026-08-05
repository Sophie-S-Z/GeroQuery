import { useEffect, useState } from "react";
import GeneView from "./components/GeneView.jsx";
import Landscape from "./components/Landscape.jsx";
import Mortality from "./components/Mortality.jsx";
import Panel from "./components/Panel.jsx";
import { json } from "./lib/db.js";
import { fmtInt } from "./lib/format.js";
import { useUrlState } from "./lib/state.js";
import { useTheme } from "./lib/theme.js";

const VIEWS = [
  ["gene", "Gene"],
  ["landscape", "Landscape"],
  ["panel", "Panel"],
  ["mortality", "Mortality"],
];

function Masthead({ meta, theme, onToggleTheme }) {
  return (
    <header className="masthead">
      <span className="wordmark">GeroQuery</span>
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
          {meta.signature_mode === "full" ? "full corpus" : "curated slice"}
        </span>
      ) : null}
      <button type="button" className="theme-toggle" onClick={onToggleTheme}>
        {theme === "dark" ? "light" : "dark"}
      </button>
    </header>
  );
}

export default function App() {
  const [state, update] = useUrlState();
  const [meta, setMeta] = useState(null);
  const [, toggleTheme, resolvedTheme] = useTheme();

  useEffect(() => {
    json("meta.json").then(setMeta).catch(() => setMeta(null));
  }, []);

  // Charts resolve CSS custom properties at draw time, so a theme change has to
  // force a redraw. Keying the view subtree on the theme is the cheapest
  // correct way to do that.
  const themeKey = resolvedTheme;

  return (
    <div className="app">
      <Masthead meta={meta} theme={resolvedTheme} onToggleTheme={toggleTheme} />

      <nav className="tabs" role="tablist" aria-label="Views">
        {VIEWS.map(([key, label]) => (
          <button
            key={key}
            type="button"
            role="tab"
            className="tab"
            aria-selected={state.view === key}
            onClick={() => update({ view: key })}
          >
            {label}
          </button>
        ))}
      </nav>

      <main key={themeKey}>
        {state.view === "gene" ? (
          <GeneView
            symbol={state.gene}
            species={state.species}
            onPick={(gene) => update({ gene })}
            onSpecies={(species) => update({ species })}
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

        {state.view === "panel" ? <Panel /> : null}
        {state.view === "mortality" ? <Mortality /> : null}
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
