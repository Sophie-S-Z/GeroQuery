import { useEffect, useState } from "react";
import { json } from "../lib/db.js";

/**
 * The corpus, auditable.
 *
 * The panel is rule-selected rather than curated, and that is the only thing
 * keeping selection bias out — so the rule is printed above the table rather
 * than described in a doc nobody opens.
 */
export default function Panel() {
  const [studies, setStudies] = useState(null);
  const [error, setError] = useState(null);

  useEffect(() => {
    let live = true;
    json("studies.json")
      .then((rows) => live && setStudies(rows))
      .catch((e) => live && setError(e.message));
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

  return (
    <>
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
          per group. Whatever passes, ships. The panel changes by changing a rule, never by
          preferring a result.
        </p>
      </div>

      <div className="panel">
        <h2 className="panel-title">
          {studies.length} contrasts · {series.size} GEO Series
        </h2>
        {shared > 0 ? (
          <p className="note" style={{ marginBottom: "0.8rem" }}>
            <strong>{shared} contrasts share a Series with another.</strong> GEO sometimes publishes
            one experiment as two DataSets — the halves of a split array — so those pairs share
            subjects. A random-effects pool assumes its inputs are independent; the Series column is
            here so that assumption can be checked rather than taken on faith.
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
