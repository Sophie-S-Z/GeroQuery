import { useEffect, useState } from "react";
import ForestPlot from "./ForestPlot.jsx";
import { contrastsForSymbol, pooledForSymbol, searchGenes } from "../lib/db.js";
import {
  VERDICT_LABEL,
  VERDICT_SUB,
  fmtEffect,
  fmtInterval,
  fmtP,
} from "../lib/format.js";

/** Debounce so a fast typist does not queue a query per keystroke. */
function useDebounced(value, delay = 160) {
  const [debounced, setDebounced] = useState(value);
  useEffect(() => {
    const timer = setTimeout(() => setDebounced(value), delay);
    return () => clearTimeout(timer);
  }, [value, delay]);
  return debounced;
}

function QueryRail({ symbol, onPick }) {
  const [term, setTerm] = useState("");
  const [results, setResults] = useState([]);
  const [error, setError] = useState(null);
  const debounced = useDebounced(term);

  useEffect(() => {
    let live = true;
    if (debounced.length < 1) {
      setResults([]);
      return undefined;
    }
    searchGenes(debounced)
      .then((rows) => live && setResults(rows))
      .catch((e) => live && setError(e.message));
    return () => {
      live = false;
    };
  }, [debounced]);

  // Collapse the per-species rows into one entry per symbol: the readout shows
  // both species side by side, so listing a symbol twice would be a false
  // distinction in the picker.
  const bySymbol = [];
  const seen = new Set();
  for (const row of results) {
    if (seen.has(row.symbol)) continue;
    seen.add(row.symbol);
    bySymbol.push(row);
  }

  return (
    <div className="rail">
      <label htmlFor="gene-search" className="panel-title">
        Gene
      </label>
      <input
        id="gene-search"
        className="search-field"
        value={term}
        placeholder="CDKN2A"
        autoComplete="off"
        spellCheck="false"
        onChange={(event) => setTerm(event.target.value)}
        onKeyDown={(event) => {
          if (event.key === "Enter" && bySymbol.length) onPick(bySymbol[0].symbol);
        }}
      />
      <p className="rail-hint">
        Symbol or alias. 41,983 genes have three or more contrasts and can be pooled.
      </p>
      {error ? <p className="error">{error}</p> : null}
      {bySymbol.length ? (
        <ul className="result-list">
          {bySymbol.map((row) => (
            <li key={row.symbol}>
              <button
                type="button"
                className="result"
                aria-current={row.symbol === symbol}
                onClick={() => onPick(row.symbol)}
              >
                <span>
                  <i className="dot" data-verdict={row.verdict} />
                  <span className="result-symbol">{row.symbol}</span>
                </span>
                <span className="result-effect">{fmtEffect(row.g)}</span>
              </button>
            </li>
          ))}
        </ul>
      ) : null}
    </div>
  );
}

export default function GeneView({ symbol, species, onPick, onSpecies }) {
  const [pooled, setPooled] = useState(null);
  const [contrasts, setContrasts] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  useEffect(() => {
    let live = true;
    setLoading(true);
    setError(null);
    Promise.all([pooledForSymbol(symbol), contrastsForSymbol(symbol, species)])
      .then(([p, c]) => {
        if (!live) return;
        setPooled(p);
        setContrasts(c);
      })
      .catch((e) => live && setError(e.message))
      .finally(() => live && setLoading(false));
    return () => {
      live = false;
    };
  }, [symbol, species]);

  const here = pooled?.find((row) => row.species === species) ?? null;
  const other = pooled?.filter((row) => row.species !== species) ?? [];
  const availableSpecies = pooled?.map((row) => row.species) ?? [];

  return (
    <div className="split">
      <QueryRail symbol={symbol} onPick={onPick} />

      <div>
        {error ? <p className="error">{error}</p> : null}

        {loading ? (
          <div className="panel" aria-busy="true">
            <div className="skeleton" style={{ width: "40%", height: "2rem" }} />
            <div className="skeleton" style={{ width: "70%" }} />
            <div className="skeleton" style={{ width: "55%" }} />
          </div>
        ) : null}

        {!loading && !here && !contrasts?.length ? (
          <div className="state">
            <p className="state-title">Not measured</p>
            <p className="mono">
              {symbol} is not in this panel, or has fewer than three contrasts in {species}.
            </p>
            <p className="note" style={{ maxWidth: "42ch", margin: "0.8rem auto 0" }}>
              That is different from finding no effect. The panel is 32 contrasts from 27 GEO
              Series; a gene absent from it was never measured here.
            </p>
          </div>
        ) : null}

        {!loading && (here || contrasts?.length) ? (
          <>
            <div className="readout-head">
              <div>
                <h1 className="gene-symbol">{symbol}</h1>
                <p className="gene-meta">
                  {species} · {contrasts?.length ?? 0} contrasts
                  {here ? ` · I² ${here.i2}% · ${here.n_tissues} tissues` : ""}
                </p>
              </div>
              {here ? (
                <p className="verdict" data-verdict={here.verdict}>
                  {VERDICT_LABEL[here.verdict]}
                  <span className="verdict-sub">{VERDICT_SUB[here.verdict]}</span>
                </p>
              ) : null}
            </div>

            {availableSpecies.length > 1 ? (
              <div className="species-switch">
                {availableSpecies.map((s) => (
                  <button key={s} type="button" aria-pressed={s === species} onClick={() => onSpecies(s)}>
                    {s}
                  </button>
                ))}
              </div>
            ) : null}

            {here ? (
              <p className="pooled-line">
                g = {fmtEffect(here.g, 3)}{" "}
                <span className="interval">{fmtInterval(here.ci_low, here.ci_high)}</span>{" "}
                <span className="interval">
                  · k = {here.k} · p = {fmtP(here.p_value)}
                </span>
              </p>
            ) : (
              <p className="note">
                Fewer than three contrasts in {species}, so no pooled estimate. The per-study
                effects below are shown without one — with one or two studies the between-study
                variance is not estimable and an interval would be a shape, not a summary.
              </p>
            )}

            <div className="panel">
              <h2 className="panel-title">Per-study evidence</h2>
              <ForestPlot contrasts={contrasts ?? []} pooled={here} symbol={symbol} />
            </div>

            {other.length ? (
              <div className="panel">
                <h2 className="panel-title">Other species</h2>
                <div className="table-scroll">
                  <table>
                    <thead>
                      <tr>
                        <th>Species</th>
                        <th className="num">g</th>
                        <th className="num">95% CI</th>
                        <th className="num">k</th>
                        <th>Verdict</th>
                      </tr>
                    </thead>
                    <tbody>
                      {other.map((row) => (
                        <tr key={row.species}>
                          <td>
                            <button
                              type="button"
                              className="link-button"
                              onClick={() => onSpecies(row.species)}
                            >
                              {row.species}
                            </button>
                          </td>
                          <td className="num">{fmtEffect(row.g)}</td>
                          <td className="num">{fmtInterval(row.ci_low, row.ci_high)}</td>
                          <td className="num">{row.k}</td>
                          <td className={row.verdict === "no_evidence" ? "null" : ""}>
                            <i className="dot" data-verdict={row.verdict} />
                            {VERDICT_LABEL[row.verdict]}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
                <p className="note" style={{ marginTop: "0.7rem" }}>
                  Species are pooled separately. A cross-species comparison is a comparison of two
                  meta-analyses, not one meta-analysis over both.
                </p>
              </div>
            ) : null}

            <div className="panel">
              <h2 className="panel-title">Contrast table</h2>
              <div className="table-scroll">
                <table>
                  <thead>
                    <tr>
                      <th>Study</th>
                      <th>Tissue</th>
                      <th>Age contrast</th>
                      <th>Sex</th>
                      <th className="num">g</th>
                      <th className="num">SE</th>
                      <th className="num">q</th>
                    </tr>
                  </thead>
                  <tbody>
                    {(contrasts ?? []).map((row) => (
                      <tr key={`${row.study_id}-${row.tissue}`}>
                        <td className="mono">{row.study_id.replace(/^GEO:/, "")}</td>
                        <td>{row.tissue ?? "—"}</td>
                        <td className="mono">{row.age_range ?? "—"}</td>
                        <td>{row.sex ?? "—"}</td>
                        <td className={`num ${row.effect_size > 0 ? "up" : "down"}`}>
                          {fmtEffect(row.effect_size, 3)}
                        </td>
                        <td className="num">{row.standard_error?.toFixed(3) ?? "—"}</td>
                        <td className="num">{fmtP(row.q_value)}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          </>
        ) : null}
      </div>
    </div>
  );
}
