import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import { css } from "../lib/theme.js";

/**
 * Weighted survival curves by dysregulation quartile, with a number-at-risk
 * table and the RMST difference called out (§4.3, §3.4).
 *
 * Hazard ratios with no survival curve is a table asking to be trusted; this is
 * the chart a reader of a mortality result expects. Four survey-weighted
 * Kaplan-Meier curves — healthiest quartile to most dysregulated — on a shared
 * time grid, with the at-risk counts beneath so the tail is readable, and the
 * restricted mean survival time (event-free years to a horizon) as the summary.
 *
 * Colour is an ordered single-hue ramp (the --seq tokens), because the quartiles
 * are ordered magnitude, not identity — never four categorical hues (§4.3).
 * Every number was computed in Python (`weighted_kaplan_meier`, `rmst`) and
 * shipped in crosslayer.json; this draws it.
 */
export default function SurvivalCurves({ survival }) {
  const host = useRef(null);
  const quartiles = survival?.quartiles ?? [];

  useEffect(() => {
    if (!host.current || quartiles.length < 2) return undefined;
    const c = css();

    // Long-form points for the step lines.
    const points = [];
    for (const q of quartiles) {
      const { grid, survival: s } = q.curve;
      for (let i = 0; i < grid.length; i++) {
        points.push({ t: grid[i], s: s[i], q: q.quartile, label: `Q${q.quartile}` });
      }
    }
    const colour = (q) => c.seq[q - 1] ?? c.ink;

    const chart = Plot.plot({
      width: 680,
      height: 380,
      marginLeft: 48,
      marginBottom: 42,
      marginRight: 34,
      style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "10px" },
      x: { label: "years of follow-up →", grid: true, nice: true },
      y: {
        label: "↑ surviving",
        domain: [0, 1],
        grid: true,
        tickFormat: (v) => `${Math.round(v * 100)}%`,
      },
      marks: [
        Plot.ruleY([0], { stroke: c.rule, strokeWidth: 1 }),
        Plot.line(points, {
          x: "t",
          y: "s",
          z: "q",
          stroke: (d) => colour(d.q),
          strokeWidth: 2,
          curve: "step-after",
        }),
        // Direct labels at the right end — four series, so label each (§5).
        Plot.text(
          quartiles.map((q) => {
            const g = q.curve.grid;
            const s = q.curve.survival;
            return { t: g[g.length - 1], s: s[s.length - 1], label: `Q${q.quartile}`, q: q.quartile };
          }),
          {
            x: "t",
            y: "s",
            text: "label",
            textAnchor: "start",
            dx: 6,
            fill: (d) => colour(d.q),
            fontSize: 10,
          },
        ),
        Plot.tip(
          points,
          Plot.pointer({
            x: "t",
            y: "s",
            title: (d) => `Q${d.q}\n${d.t.toFixed(1)} y · ${(d.s * 100).toFixed(1)}% surviving`,
          }),
        ),
      ],
    });
    host.current.replaceChildren(chart);
    return () => chart.remove();
  }, [quartiles]);

  if (quartiles.length < 2) return null;

  const horizons = survival.horizons ?? [10, 15];
  const byQ = Object.fromEntries(quartiles.map((q) => [q.quartile, q]));
  const healthiest = byQ[1];
  const worst = byQ[4] ?? quartiles[quartiles.length - 1];
  const ticks = survival.grid.filter((t, i) => i % Math.ceil(survival.grid.length / 6) === 0);
  const nearestIdx = (t) => {
    const g = quartiles[0].curve.grid;
    let best = 0;
    for (let i = 1; i < g.length; i++) if (Math.abs(g[i] - t) < Math.abs(g[best] - t)) best = i;
    return best;
  };

  return (
    <div className="panel">
      <h2 className="panel-title">Survival by dysregulation quartile</h2>
      <div className="chart" ref={host} role="img" aria-label="Weighted survival curves by dysregulation quartile" />

      <div className="table-scroll" style={{ marginTop: "0.6rem" }}>
        <table className="at-risk">
          <thead>
            <tr>
              <th>At risk</th>
              {ticks.map((t) => (
                <th key={t} className="num">
                  {Math.round(t)}y
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {quartiles.map((q) => (
              <tr key={q.quartile}>
                <td>
                  <i
                    className="swatch"
                    style={{ background: `var(--seq-${q.quartile})`, width: "9px", height: "9px" }}
                  />{" "}
                  Q{q.quartile}
                </td>
                {ticks.map((t) => (
                  <td key={t} className="num">
                    {q.curve.n_at_risk[nearestIdx(t)]}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {healthiest && worst ? (
        <p className="note" style={{ marginTop: "0.8rem", maxWidth: "70ch" }}>
          {horizons.map((h) => {
            const diff = healthiest.rmst[String(h)] - worst.rmst[String(h)];
            return (
              <span key={h} style={{ display: "block" }}>
                By <strong className="mono">{h}</strong> years, the most-dysregulated quartile has
                lived <strong className="mono">{diff.toFixed(2)}</strong> fewer event-free years
                than the healthiest ({worst.rmst[String(h)].toFixed(2)} vs{" "}
                {healthiest.rmst[String(h)].toFixed(2)}) — restricted mean survival time, which a
                hazard ratio cannot say.
              </span>
            );
          })}
        </p>
      ) : null}
    </div>
  );
}
