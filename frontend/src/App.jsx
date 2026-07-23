import React, { useEffect, useState } from "react";
import { api } from "./api.js";

const EXAMPLES = ["CDKN2A", "LMNB1", "KLOTHO", "GDF15", "SIRT1", "FOXO3"];

function EffectBar({ effect, max }) {
  const pct = Math.min(100, (Math.abs(effect) / max) * 100);
  const up = effect >= 0;
  return (
    <div className="effect-track">
      <div className="effect-mid" />
      <div
        className={`effect-fill ${up ? "up" : "down"}`}
        style={{ width: `${pct / 2}%`, [up ? "left" : "right"]: "50%" }}
      />
    </div>
  );
}

export default function App() {
  const [query, setQuery] = useState("CDKN2A");
  const [species, setSpecies] = useState("");
  const [card, setCard] = useState(null);
  const [status, setStatus] = useState("idle");
  const [error, setError] = useState(null);
  const [version, setVersion] = useState(null);

  useEffect(() => {
    api.version().then(setVersion).catch(() => {});
  }, []);

  async function search(q = query, sp = species) {
    setStatus("loading");
    setError(null);
    try {
      const data = await api.geneCard(q, sp);
      setCard(data);
      setStatus("done");
    } catch (e) {
      setError(e.message);
      setStatus("error");
    }
  }

  useEffect(() => {
    search();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const metas = card?.meta_signatures ?? [];
  const maxEffect = Math.max(0.5, ...metas.map((m) => Math.abs(m.pooled_effect)));

  return (
    <div className="page">
      <header>
        <div className="brand">
          <span className="logo">🧬</span>
          <div>
            <h1>GeroQuery</h1>
            <p className="sub">
              Multi-omic · cross-species · clock- &amp; resilience-aware aging data
            </p>
          </div>
        </div>
        {version && (
          <span className="version">data {version.data_version}</span>
        )}
      </header>

      <form
        className="search"
        onSubmit={(e) => {
          e.preventDefault();
          search();
        }}
      >
        <input
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="Gene symbol, Ensembl, Entrez, or alias…"
          aria-label="Gene query"
        />
        <select value={species} onChange={(e) => setSpecies(e.target.value)} aria-label="Species">
          <option value="">both species</option>
          <option value="human">human</option>
          <option value="mouse">mouse</option>
        </select>
        <button type="submit">Search</button>
      </form>

      <div className="examples">
        {EXAMPLES.map((g) => (
          <button key={g} className="chip" onClick={() => { setQuery(g); search(g, species); }}>
            {g}
          </button>
        ))}
      </div>

      {status === "loading" && <div className="notice">Loading…</div>}
      {status === "error" && <div className="notice error">{error}</div>}

      {status === "done" && card && (
        <main>
          <section className="card gene-head">
            <h2>
              {card.gene.symbol} <span className="muted">{card.gene.name}</span>
            </h2>
            <div className="ids">
              <code>{card.gene.canonical_id}</code>
              <span>Entrez {card.gene.entrez ?? "—"}</span>
              <span>UniProt {card.gene.uniprot ?? "—"}</span>
              <span>{card.gene.species}</span>
            </div>
          </section>

          <section className="card">
            <h3>Pooled aging effect</h3>
            {metas.length === 0 ? (
              <p className="muted">No harmonized signatures in the bundled slice.</p>
            ) : (
              <table>
                <thead>
                  <tr>
                    <th>omic</th>
                    <th>species</th>
                    <th>effect (g)</th>
                    <th></th>
                    <th>I²</th>
                    <th>studies</th>
                  </tr>
                </thead>
                <tbody>
                  {metas.map((m, i) => (
                    <tr key={i}>
                      <td>{m.omic_layer}</td>
                      <td>{m.species}</td>
                      <td className={m.direction === "up" ? "up" : "down"}>
                        {m.pooled_effect > 0 ? "+" : ""}
                        {m.pooled_effect}
                      </td>
                      <td className="barcell">
                        <EffectBar effect={m.pooled_effect} max={maxEffect} />
                      </td>
                      <td>{m.heterogeneity_i2}%</td>
                      <td>{m.n_studies}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            )}
          </section>

          <div className="grid2">
            <section className="card">
              <h3>Curated knowledge</h3>
              {card.curated_flags.length === 0 ? (
                <p className="muted">None in slice.</p>
              ) : (
                <ul className="flags">
                  {card.curated_flags.map((f, i) => (
                    <li key={i}>
                      <span className="db">{f.database}</span> {f.assertion}
                    </li>
                  ))}
                </ul>
              )}
            </section>

            <section className="card">
              <h3>Linked interventions</h3>
              {card.interventions.length === 0 ? (
                <p className="muted">None linked in slice.</p>
              ) : (
                <ul className="flags">
                  {card.interventions.map((iv, i) => (
                    <li key={i}>
                      <span className="db">{iv.source}</span> {iv.name}
                      {iv.lifespan_effect_pct != null && (
                        <span className="eff"> +{iv.lifespan_effect_pct}% lifespan</span>
                      )}
                    </li>
                  ))}
                </ul>
              )}
            </section>
          </div>
        </main>
      )}

      <footer>
        Calls the same <code>/v1</code> API as every other client. Open-source · MIT.
      </footer>
    </div>
  );
}
