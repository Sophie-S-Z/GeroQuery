import React, { useEffect, useState } from "react";
import { api, REAL_CLINICAL_DATASET } from "./api.js";

const EXAMPLES = ["CDKN2A", "LMNB1", "KLOTHO", "GDF15", "SIRT1", "FOXO3"];

const fmt = (x, digits = 5) =>
  x == null || Number.isNaN(x) ? "—" : `${x > 0 ? "+" : ""}${x.toFixed(digits)}`;

/** One critical-slowing-down indicator, reported as evidence rather than sign.
 *  A positive slope alone is not a finding — the interval has to exclude zero. */
function TrendRow({ label, evidence }) {
  if (!evidence) return null;
  const { slope, ci_low, ci_high, supported } = evidence;
  return (
    <tr>
      <td>{label}</td>
      <td className={supported ? "up" : "muted"}>{fmt(slope)}</td>
      <td className="muted">
        [{fmt(ci_low)}, {fmt(ci_high)}]
      </td>
      <td>
        <span className={`badge ${supported ? "badge-yes" : "badge-no"}`}>
          {supported ? "supported" : "no evidence"}
        </span>
      </td>
    </tr>
  );
}

function ResiliencePanel({ csd, error, datasetNote }) {
  if (error) {
    return (
      <section className="card">
        <h3>Resilience — real NHANES</h3>
        <p className="notice error">{error}</p>
      </section>
    );
  }
  if (!csd) {
    return (
      <section className="card">
        <h3>Resilience — real NHANES</h3>
        <p className="muted">Loading…</p>
      </section>
    );
  }
  return (
    <section className="card">
      <h3>
        Resilience — real NHANES <span className="tag tag-real">real data</span>
      </h3>
      <p className="muted small">
        Critical slowing down across {csd.strata_midpoints?.length ?? 0} age strata,
        n&nbsp;=&nbsp;{csd.n_samples} subjects. Trends are gated on a subject-level
        bootstrap CI, not on the sign of the slope.
      </p>
      <table>
        <thead>
          <tr>
            <th>indicator</th>
            <th>slope / yr</th>
            <th>95% CI</th>
            <th>verdict</th>
          </tr>
        </thead>
        <tbody>
          <TrendRow label="health-state variance" evidence={csd.variance_evidence} />
          <TrendRow label="marker cross-correlation" evidence={csd.crosscorr_evidence} />
        </tbody>
      </table>
      <p className={`conclusion ${csd.resilience_declines ? "up" : "down"}`}>
        {csd.resilience_declines
          ? "Both early-warning signals are supported."
          : "Only one of the two early-warning signals is supported — this is reported as a partial result, not as resilience decline."}
      </p>
      {datasetNote && <p className="muted small">{datasetNote}</p>}
    </section>
  );
}

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
  const [csd, setCsd] = useState(null);
  const [csdError, setCsdError] = useState(null);
  const [datasetNote, setDatasetNote] = useState(null);

  useEffect(() => {
    api.version().then(setVersion).catch(() => {});
    api
      .csd(REAL_CLINICAL_DATASET)
      .then(setCsd)
      .catch((e) => setCsdError(e.message));
    // The store records whether the full cohort or the offline sample was built.
    // Surfacing it stops a 600-row sample estimate being read as the headline.
    api
      .datasets()
      .then(({ datasets }) => {
        const row = datasets?.find((d) => d.dataset_id === REAL_CLINICAL_DATASET);
        if (row) setDatasetNote(row.description);
      })
      .catch(() => {});
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

          <ResiliencePanel csd={csd} error={csdError} datasetNote={datasetNote} />

          <section className="card">
            <h3>
              Pooled aging effect{" "}
              <span className="tag tag-real">real data</span>
            </h3>
            <p className="muted small">
              Random-effects pool (DerSimonian&ndash;Laird) over young-vs-old contrasts
              from checksum-pinned NCBI GEO DataSets. Read the interval, not the point:
              with 3&ndash;15 samples per study, a confidence interval spanning zero
              means <em>not detected by this panel</em>, not <em>absent</em>.
            </p>
            {metas.length === 0 ? (
              <p className="muted">No GEO contrast in this panel covers this gene.</p>
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
