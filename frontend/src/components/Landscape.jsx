import { useEffect, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import { css } from "../lib/theme.js";
import { effectDistribution, panelSummary, topGenes } from "../lib/db.js";
import { fmtEffect, fmtInterval, fmtLfsr } from "../lib/format.js";

/**
 * The whole corpus at once.
 *
 * Until now there was no way to *see* 485,905 rows, only to query into them.
 * This is that view: every pooled estimate binned by effect size and stacked by
 * verdict, so the shape of the evidence is visible — including how much of it
 * straddles zero, which is most of it.
 */
export default function Landscape({ species, direction, onPick, onSpecies, onDirection }) {
  const host = useRef(null);
  const [summary, setSummary] = useState(null);
  const [rows, setRows] = useState([]);
  const [bins, setBins] = useState([]);
  const [error, setError] = useState(null);

  useEffect(() => {
    let live = true;
    setError(null);
    Promise.all([
      panelSummary(),
      effectDistribution(species),
      topGenes({ species, direction }),
    ])
      .then(([s, b, t]) => {
        if (!live) return;
        setSummary(s);
        setBins(b);
        setRows(t);
      })
      .catch((e) => live && setError(e.message));
    return () => {
      live = false;
    };
  }, [species, direction]);

  useEffect(() => {
    if (!host.current || !bins.length) return undefined;
    const c = css();
    const chart = Plot.plot({
      width: 860,
      height: 260,
      marginLeft: 54,
      marginBottom: 38,
      style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "10px" },
      x: { label: "Pooled effect (Hedges' g) →", nice: true },
      y: { label: "↑ genes", grid: true },
      color: {
        domain: ["decreases", "no_evidence", "increases"],
        range: [c.down, c.null_, c.up],
      },
      marks: [
        Plot.rectY(bins, {
          x: "bin",
          interval: 0.1,
          y: "n",
          fill: "verdict",
          title: (d) => `g ≈ ${d.bin.toFixed(1)}\n${d.n} genes\n${d.verdict}`,
        }),
        // Zero again, and again heavier than the grid.
        Plot.ruleX([0], { stroke: c.ink, strokeWidth: 1.75 }),
        Plot.ruleY([0], { stroke: c.rule }),
      ],
    });
    host.current.replaceChildren(chart);
    return () => chart.remove();
  }, [bins]);

  const here = summary?.find((row) => row.species === species);

  return (
    <>
      <div className="panel">
        <h2 className="panel-title">Distribution of pooled effects — {species}</h2>
        {error ? <p className="error">{error}</p> : null}
        <div className="species-switch">
          {["human", "mouse"].map((s) => (
            <button
              key={s}
              type="button"
              aria-pressed={species === s}
              onClick={() => onSpecies(s)}
            >
              {s}
            </button>
          ))}
        </div>
        <div className="chart" ref={host} role="img" aria-label="Effect size distribution" />
        {here ? (
          <p className="note" style={{ marginTop: "0.8rem" }}>
            <strong className="mono">{here.n_pooled.toLocaleString()}</strong> genes pooled from
            three or more contrasts.{" "}
            <span className="up mono">{here.n_up.toLocaleString()}</span> rise with age,{" "}
            <span className="down mono">{here.n_down.toLocaleString()}</span> fall, and{" "}
            <span className="null mono">{here.n_null.toLocaleString()}</span> have an interval that
            crosses zero. Median interval width is{" "}
            <strong className="mono">{here.median_ci_width}</strong> — wide enough that this panel
            detects large consistent effects and cannot rule out moderate ones.
          </p>
        ) : null}
      </div>

      <div className="panel">
        <h2 className="panel-title">
          Strongest replicated effects — {direction === "up" ? "rising" : "falling"} with age
        </h2>
        <div className="species-switch">
          {[
            ["up", "rising"],
            ["down", "falling"],
          ].map(([key, label]) => (
            <button
              key={key}
              type="button"
              aria-pressed={direction === key}
              onClick={() => onDirection(key)}
            >
              {label}
            </button>
          ))}
        </div>
        <div className="table-scroll">
          <table>
            <thead>
              <tr>
                <th>Gene</th>
                <th className="num">g</th>
                <th className="num">shrunk</th>
                <th className="num">95% CI</th>
                <th className="num">lfsr</th>
                <th className="num">k</th>
                <th className="num">I²</th>
                <th className="num">tissues</th>
              </tr>
            </thead>
            <tbody>
              {rows.map((row) => (
                <tr key={row.gene_id}>
                  <td>
                    <button
                      type="button"
                      className="link-button mono"
                      onClick={() => onPick(row.symbol)}
                    >
                      {row.symbol}
                    </button>
                  </td>
                  <td className={`num ${direction === "up" ? "up" : "down"}`}>
                    {fmtEffect(row.g)}
                  </td>
                  <td className="num">{fmtEffect(row.g_shrunk)}</td>
                  <td className="num">{fmtInterval(row.ci_low, row.ci_high)}</td>
                  {/* Faded above 0.2: a gene whose direction is more than
                      one-in-five likely to be wrong should not read with the
                      same weight as one that is settled. */}
                  <td className="num" style={row.lfsr > 0.2 ? { color: "var(--ink-faint)" } : null}>
                    {fmtLfsr(row.lfsr)}
                  </td>
                  <td className="num">{row.k}</td>
                  <td className="num">{row.i2}%</td>
                  <td className="num">{row.n_tissues}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <p className="note" style={{ marginTop: "0.7rem" }}>
          Ranked by the <strong>shrunken</strong> effect among genes whose interval excludes zero,
          not by the raw one. Ranking by raw <span className="mono">g</span> selects on
          significance and then sorts on an extreme of a noisy statistic, which reliably puts the
          genes most likely to shrink on replication at the top: under the old ordering the mouse
          table led with <span className="mono">Kif21a</span>, whose interval excludes zero and
          whose posterior gives its direction a <strong>46% chance of being wrong</strong>. The{" "}
          <span className="mono">lfsr</span> column is that probability, and it is faded where it
          exceeds 0.2. Rank order is still not a claim about importance — a large effect measured
          in one tissue is not stronger evidence than a moderate one measured in six.
        </p>
      </div>
    </>
  );
}
