import { useEffect, useMemo, useState } from "react";
import CorpusInfluence from "./CorpusInfluence.jsx";
import { json } from "../lib/db.js";

/**
 * The corpus, auditable.
 *
 * The panel is rule-selected rather than curated, and that is the only thing
 * keeping selection bias out — so the rule is printed above the table rather
 * than described in a doc nobody opens. Above it now sits the corpus-influence
 * chart (§4.5): the same contrasts, ranked by how much each one moves the whole
 * published corpus — provenance made measurable.
 */
export default function Panel() {
  const [studies, setStudies] = useState(null);
  const [influence, setInfluence] = useState(null);
  const [error, setError] = useState(null);
  // Memoised above the early returns (hook order) so the chart, which depends
  // on it, redraws when the data changes rather than on every render.
  const byStudy = useMemo(
    () => Object.fromEntries((studies ?? []).map((s) => [s.study_id, s])),
    [studies],
  );

  useEffect(() => {
    let live = true;
    json("studies.json")
      .then((rows) => live && setStudies(rows))
      .catch((e) => live && setError(e.message));
    // Optional: a clean checkout that has not run the influence enrichment
    // simply shows the panel without the ranking, rather than failing.
    json("influence_corpus.json")
      .then((rows) => live && setInfluence(rows))
      .catch(() => live && setInfluence([]));
    return () => {
      live = false;
    };
  }, []);

  if (error) return <p className="error">{error}</p>;
  if (!studies) {
    return (
      <div className="panel">
        <div className="skeleton" style={{ width: "60%" }} />
        <div className="skeleton" style={{ width: "45%" }} />
      </div>
    );
  }

  // Count only the contrasts that actually declare a Series. `series_id` is
  // nullable (it renders as "—" in the table below), and subtracting the
  // distinct-Series count from *all* contrasts charged every null one as
  // sharing with another — inventing a non-independence warning on the page
  // whose whole job is letting a reader audit that assumption. Today's panel
  // happens to have no nulls, so the old arithmetic was right by luck; the
  // panel is rebuilt monthly from upstream GEO.
  const withSeries = studies.filter((s) => s.series_id);
  const series = new Set(withSeries.map((s) => s.series_id));
  const shared = withSeries.length - series.size;

  const influenceByStudy = Object.fromEntries((influence ?? []).map((r) => [r.study_id, r]));
  const hasInfluence = influence && influence.length > 0;

  return (
    <>
      {hasInfluence ? (
        <div className="panel">
          <h2 className="panel-title">Which contrast bends the corpus</h2>
          <p className="note" style={{ marginBottom: "0.8rem" }}>
            Each bar is how far the pooled estimate moves, on average, when that one contrast is
            dropped from a gene it measures. A contrast covers only the genes its platform measures —
            between {Math.min(...influence.map((r) => r.n_genes)).toLocaleString()} and{" "}
            {Math.max(...influence.map((r) => r.n_genes)).toLocaleString()} here — so the average,
            not the total, is shown: a total would mostly rank platform coverage. The species are drawn apart
            because a mouse gene is pooled from fewer, noisier contrasts, so dropping one moves it
            further by construction — compare contrasts within a species. This is a question about
            the panel itself, and nothing in a curated database can answer it.
          </p>
          <CorpusInfluence ranking={influence} byStudy={byStudy} />
        </div>
      ) : null}

      <div className="panel">
        <h2 className="panel-title">How this panel was selected</h2>
        <p className="mono" style={{ fontSize: "0.76rem", lineHeight: 1.65, margin: "0 0 0.8rem" }}>
          age[Subset Variable Type] AND gds[Filter]
          <br />
          AND (Homo sapiens[Organism] OR Mus musculus[Organism])
        </p>
        <p className="note">
          Every GEO DataSet declaring an age subset variable, then an adult young-versus-old split
          restricted to the control arm of every other subset variable, with at least three samples
          per group.
        </p>
        <p className="note">
          Then every bulk RNA-seq sample in ARCHS4 (one uniform alignment of all public GEO RNA-seq)
          whose characteristics record an age in the same young or old band. Within a series,
          samples are grouped by tissue, cell type and every other recorded attribute, so the two
          arms differ only in age; treatments, genotypes, diets and diseases must sit at their
          control level, any unrecognised attribute excludes the sample, and a series that profiles
          any tumour is left out whole. At least three samples per group.
        </p>
        <p className="note">
          Whatever passes, ships. The panel changes by changing a rule, never by preferring a
          result.
        </p>
      </div>

      <div className="panel">
        <h2 className="panel-title">
          {studies.length} contrasts · {series.size} GEO Series
        </h2>
        {shared > 0 ? (
          <p className="note" style={{ marginBottom: "0.8rem" }}>
            <strong>{shared} contrasts share a Series with another.</strong> One experiment can yield
            several contrasts — the two halves of a split array, or several tissues and cell types
            from the same animals — and those share subjects. Each Series is therefore combined
            into one estimate before pooling, so a gene&rsquo;s k counts independent Series, not
            contrasts; the Series column shows which contrasts were combined.
          </p>
        ) : null}
        <div className="table-scroll">
          <table>
            <thead>
              <tr>
                <th>Study</th>
                <th>Series</th>
                <th>Species</th>
                <th>Tissue</th>
                <th className="num">n</th>
                {hasInfluence ? <th className="num">mean shift</th> : null}
                <th>Platform</th>
                <th>Source</th>
              </tr>
            </thead>
            <tbody>
              {studies.map((study) => (
                <tr key={study.study_id}>
                  <td className="mono">
                    {study.url ? (
                      <a href={study.url} target="_blank" rel="noreferrer">
                        {study.study_id.replace(/^GEO:/, "")}
                      </a>
                    ) : (
                      study.study_id.replace(/^GEO:/, "")
                    )}
                  </td>
                  <td className="mono">{study.series_id ?? "—"}</td>
                  <td>{study.species}</td>
                  <td>{study.tissue ?? "—"}</td>
                  <td className="num">{study.sample_size ?? "—"}</td>
                  {hasInfluence ? (
                    <td className="num" title="mean leave-one-out shift per gene this contrast measures">
                      {influenceByStudy[study.study_id]
                        ? influenceByStudy[study.study_id].mean_abs_shift.toFixed(3)
                        : "—"}
                    </td>
                  ) : null}
                  <td className="mono">{study.processing_method ?? "—"}</td>
                  <td className="note">{study.source}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </>
  );
}
