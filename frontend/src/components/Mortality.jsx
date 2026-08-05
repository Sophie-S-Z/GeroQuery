import { useEffect, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import { css } from "../lib/theme.js";
import { json } from "../lib/db.js";
import { fmtHR, fmtHRInterval, fmtInt, fmtP } from "../lib/format.js";

/**
 * The cross-layer result: clocks, health state, and death on the same people.
 *
 * The chart's reference line is 1.0 rather than 0 — hazard ratios are
 * multiplicative — but it does the same job as the zero rule elsewhere, and is
 * drawn the same weight for the same reason.
 */
export default function Mortality() {
  const host = useRef(null);
  const [data, setData] = useState(null);
  const [error, setError] = useState(null);

  useEffect(() => {
    let live = true;
    json("crosslayer.json")
      .then((payload) => live && setData(payload))
      .catch((e) => live && setError(e.message));
    return () => {
      live = false;
    };
  }, []);

  useEffect(() => {
    if (!host.current || !data?.available) return undefined;
    const c = css();
    const rows = data.clocks.map((row) => ({
      ...row,
      // Strip only the "Mort" suffix. An earlier version also stripped "Age",
      // via /Age$|Mort$|AgeMort$/ — but regex alternation picks the leftmost
      // position first, so "GrimAgeMort" matched `AgeMort$` at index 4 and
      // became "Grim" while "GrimAge2Mort" only matched `Mort$` and became
      // "GrimAge2". Two rows of the same clock family, labelled inconsistently.
      label: row.clock.replace(/Mort$/, ""),
    }));

    const chart = Plot.plot({
      width: 760,
      height: rows.length * 30 + 80,
      marginLeft: 108,
      marginBottom: 42,
      style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "10px" },
      x: { label: "Hazard ratio per SD of age acceleration →", nice: true, grid: true },
      y: { label: null, domain: rows.map((r) => r.label) },
      marks: [
        Plot.ruleX([1], { stroke: c.ink, strokeWidth: 1.75 }),
        Plot.ruleY(rows, {
          y: "label",
          x1: "ci_low",
          x2: "ci_high",
          stroke: (d) => (d.excludes_null ? c.up : c.null_),
          strokeWidth: 1.5,
          strokeLinecap: "round",
        }),
        Plot.dot(rows, {
          y: "label",
          x: "hr",
          fill: (d) => (d.excludes_null ? c.up : c.paper),
          stroke: (d) => (d.excludes_null ? c.up : c.null_),
          strokeWidth: 1.5,
          r: 3.8,
          title: (d) =>
            `${d.clock}\ntrained on ${d.predicts.replace(/_/g, " ")}\n` +
            `HR ${fmtHR(d.hr)} ${fmtHRInterval(d.ci_low, d.ci_high)}\n${d.note}`,
        }),
      ],
    });
    host.current.replaceChildren(chart);
    return () => chart.remove();
  }, [data]);

  if (error) return <p className="error">{error}</p>;
  if (!data) return <div className="panel"><div className="skeleton" style={{ width: "50%" }} /></div>;

  if (!data.available) {
    return (
      <div className="state">
        <p className="state-title">Cohort not built</p>
        <p className="mono">Run `make crosslayer` and rebuild the site data.</p>
        <p className="note" style={{ maxWidth: "46ch", margin: "0.8rem auto 0" }}>
          The cross-layer table is reproducible from the pinned manifest, so it is not committed.
          {data.reason ? ` ${data.reason}` : ""}
        </p>
      </div>
    );
  }

  const best = data.clocks[0];

  return (
    <>
      <div className="panel">
        <h2 className="panel-title">
          NHANES 1999–2002 · clocks, health state, and mortality on the same people
        </h2>
        <p className="note" style={{ marginBottom: "1rem" }}>
          <strong className="mono">{fmtInt(data.n_subjects)}</strong> adults aged{" "}
          {data.age_min}–{data.age_max}, <strong className="mono">{fmtInt(data.n_deaths)}</strong>{" "}
          deaths over a median{" "}
          <strong className="mono">{data.followup_years_median}</strong> years. Twelve DNA
          methylation clocks computed by NCHS, joined on <span className="mono">SEQN</span> to the
          six-marker health state and to the 2019 linked mortality file.
        </p>
        <div className="chart" ref={host} role="img" aria-label="Hazard ratios by clock" />
        <div className="legend">
          <span>
            <i className="swatch" style={{ background: "var(--up)" }} /> interval excludes 1.0
          </span>
          <span>
            <i className="swatch" style={{ background: "var(--null)" }} /> interval crosses 1.0
          </span>
        </div>
      </div>

      <div className="panel">
        <h2 className="panel-title">Does either add anything over the other?</h2>
        <div className="table-scroll">
          <table>
            <thead>
              <tr>
                <th>Clock</th>
                <th>Trained on</th>
                <th className="num">HR / SD</th>
                <th className="num">95% CI</th>
                <th className="num">C clock</th>
                <th className="num">C joint</th>
                <th className="num">health state adds</th>
                <th className="num">clock adds</th>
              </tr>
            </thead>
            <tbody>
              {data.clocks.map((row) => (
                <tr key={row.clock}>
                  <td className="mono">{row.clock}</td>
                  <td className="note">{row.predicts.replace(/_/g, " ")}</td>
                  <td className={`num ${row.excludes_null ? "up" : "null"}`}>{fmtHR(row.hr)}</td>
                  <td className="num">{fmtHRInterval(row.ci_low, row.ci_high)}</td>
                  <td className="num">{row.c_clock}</td>
                  <td className="num">{row.c_joint}</td>
                  <td className="num">{fmtP(row.dysregulation_adds_p)}</td>
                  <td className={`num ${row.clock_adds_p > 0.05 ? "null" : ""}`}>
                    {fmtP(row.clock_adds_p)}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <p className="note" style={{ marginTop: "0.9rem" }}>
          The last column is the one worth reading twice. For{" "}
          <span className="mono">SkinBloodAge</span>, <span className="mono">LinAge</span> and{" "}
          <span className="mono">WeidnerAge</span> it is above 0.05: there is no evidence those
          clocks tell you anything six routine blood tests do not.
        </p>
      </div>

      <div className="panel">
        <h2 className="panel-title">The number that deflates the p-values</h2>
        <p className="note">
          Age and sex alone reach C ={" "}
          <strong className="mono">{best?.c_baseline}</strong>. The health state alone reaches{" "}
          <strong className="mono">{best?.c_dysregulation}</strong> — better than nine of the ten
          clocks. Together with the strongest clock, C ={" "}
          <strong className="mono">{best?.c_joint}</strong>.
        </p>
        <p className="note" style={{ marginTop: "0.6rem" }}>
          So all of this biology adds{" "}
          <strong className="mono">
            {best ? (best.c_joint - best.c_baseline).toFixed(3) : "—"}
          </strong>{" "}
          of concordance over knowing someone&rsquo;s age and sex. With{" "}
          {fmtInt(data.n_deaths)} deaths a likelihood-ratio test detects effects far below the
          size that would change a decision about a person. Read the interval, not the exponent.
        </p>
      </div>

      <div className="panel">
        <h2 className="panel-title">Limitations</h2>
        <ul className="caveats">
          {data.caveats.map((caveat) => (
            <li key={caveat}>{caveat}</li>
          ))}
          <li>
            DNAm-predicted sex disagrees with reported sex on{" "}
            <span className="mono">{data.sex_qc.n_discordant}</span> of{" "}
            <span className="mono">{data.sex_qc.n_compared}</span> samples (
            {(data.sex_qc.rate * 100).toFixed(1)}%). Reported rather than dropped: a rise in that
            rate is the cheapest signal that samples were mixed up upstream.
          </li>
        </ul>
      </div>
    </>
  );
}
