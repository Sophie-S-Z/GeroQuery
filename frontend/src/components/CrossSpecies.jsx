import { useEffect, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import { json } from "../lib/db.js";
import { css } from "../lib/theme.js";
import { fmtEffect } from "../lib/format.js";

/**
 * Cross-species scatter (§4.2) — the view most likely to produce a finding.
 *
 * Human pooled effect against mouse pooled effect, one point per ortholog group
 * present in both species. The primary mark is the **y = x identity line**: the
 * question is agreement, not slope. A Deming/orthogonal fit is drawn too —
 * computed in Python (`deming_regression`) and shipped, because both axes carry
 * measurement error and an OLS line would be biased flat and understate
 * conservation. Points are coloured by whether the two species agree, using the
 * polarity palette; concordance is read off the verdicts, not the raw sign.
 *
 * The honest headline here is how *few* pairs resolve in both species at once —
 * that is the finding, and the chart shows it rather than hiding it.
 */
export default function CrossSpecies({ onPick }) {
  const host = useRef(null);
  const [data, setData] = useState(null);
  const [error, setError] = useState(null);

  useEffect(() => {
    let live = true;
    json("crossspecies.json")
      .then((d) => live && setData(d))
      .catch((e) => live && setError(e.message));
    return () => {
      live = false;
    };
  }, []);

  const quadrant = (p) => {
    const { verdict_human: h, verdict_mouse: m } = p;
    if (h === "no_evidence" || m === "no_evidence") return "unresolved";
    if (h === m) return h === "increases" ? "both_up" : "both_down";
    return "discordant";
  };

  useEffect(() => {
    if (!host.current || !data?.available) return undefined;
    const c = css();
    const colour = {
      both_up: c.up,
      both_down: c.down,
      discordant: c.ink,
      unresolved: c.null_,
    };

    const pts = data.pairs.map((p) => ({ ...p, q: quadrant(p) }));
    const xs = pts.map((p) => p.g_human);
    const ys = pts.map((p) => p.g_mouse);
    const lo = Math.min(-1, ...xs, ...ys) - 0.1;
    const hi = Math.max(1, ...xs, ...ys) + 0.1;
    const { slope, intercept } = data.deming;

    // Draw the resolved points last so they sit above the grey cloud.
    const order = { unresolved: 0, both_up: 1, both_down: 1, discordant: 2 };
    pts.sort((a, b) => order[a.q] - order[b.q]);

    const chart = Plot.plot({
      width: 680,
      height: 680,
      marginLeft: 54,
      marginBottom: 48,
      marginTop: 12,
      marginRight: 16,
      style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "11px" },
      x: { domain: [lo, hi], label: "human pooled g →", grid: true, nice: true },
      y: { domain: [lo, hi], label: "↑ mouse pooled g", grid: true, nice: true },
      marks: [
        Plot.ruleX([0], { stroke: c.rule, strokeWidth: 1 }),
        Plot.ruleY([0], { stroke: c.rule, strokeWidth: 1 }),
        // Identity: the reference the question is actually about.
        Plot.link([{ x1: lo, y1: lo, x2: hi, y2: hi }], {
          x1: "x1",
          y1: "y1",
          x2: "x2",
          y2: "y2",
          stroke: c.ink,
          strokeWidth: 1.5,
        }),
        // Deming fit, dashed, drawn from Python's slope/intercept.
        Plot.link([{ x1: lo, y1: slope * lo + intercept, x2: hi, y2: slope * hi + intercept }], {
          x1: "x1",
          y1: "y1",
          x2: "x2",
          y2: "y2",
          stroke: c.inkFaint,
          strokeWidth: 1.5,
          strokeDasharray: "5,4",
        }),
        Plot.dot(pts, {
          x: "g_human",
          y: "g_mouse",
          r: (d) => (d.q === "unresolved" ? 2.3 : 4),
          fill: (d) => colour[d.q],
          fillOpacity: (d) => (d.q === "unresolved" ? 0.45 : 0.95),
          stroke: (d) => (d.q === "unresolved" ? "none" : c.paper),
          strokeWidth: 0.8,
          title: (d) =>
            `${d.symbol_human} / ${d.symbol_mouse}\n` +
            `human g = ${d.g_human.toFixed(3)} (${d.verdict_human.replace("_", " ")})\n` +
            `mouse g = ${d.g_mouse.toFixed(3)} (${d.verdict_mouse.replace("_", " ")})`,
        }),
      ],
    });

    host.current.replaceChildren(chart);
    return () => chart.remove();
  }, [data]);

  if (error) return <p className="error">{error}</p>;
  if (!data) {
    return (
      <div className="panel">
        <div className="skeleton" style={{ width: "50%" }} />
        <div className="skeleton" style={{ width: "40%" }} />
      </div>
    );
  }
  if (!data.available) {
    return (
      <div className="state">
        <p className="state-title">Not built</p>
        <p className="note">
          The cross-species pairing needs the exporter&rsquo;s ortholog-group column; this build
          does not carry it.
        </p>
      </div>
    );
  }

  const q = data.quadrants;
  const resolved = data.pairs.filter(
    (p) => p.verdict_human !== "no_evidence" && p.verdict_mouse !== "no_evidence",
  );
  const d = data.deming;

  return (
    <div>
      <div className="readout-head">
        <div>
          <h1 className="section-title">Human vs mouse</h1>
          <p className="gene-meta">
            {data.n_pairs.toLocaleString()} ortholog groups measured in both species
          </p>
        </div>
      </div>

      <p className="note" style={{ maxWidth: "70ch" }}>
        Each point is one gene&rsquo;s pooled age effect in human against its ortholog&rsquo;s in
        mouse. The solid line is <strong>y = x</strong> — perfect agreement — and the question is how
        close the cloud sits to it, not its slope. The dashed line is a Deming (orthogonal) fit,
        computed in Python because both axes carry error and an ordinary least-squares line would be
        biased toward flat and understate conservation.
      </p>

      <div className="panel">
        <div className="chart" ref={host} role="img" aria-label="Human versus mouse pooled effects" />
        <div className="legend">
          <span>
            <i className="swatch" style={{ background: "var(--up)", height: "9px", width: "9px", borderRadius: "50%" }} />{" "}
            both rise ({q.both_up})
          </span>
          <span>
            <i className="swatch" style={{ background: "var(--down)", height: "9px", width: "9px", borderRadius: "50%" }} />{" "}
            both fall ({q.both_down})
          </span>
          <span>
            <i className="swatch" style={{ background: "var(--ink)", height: "9px", width: "9px", borderRadius: "50%" }} />{" "}
            discordant ({q.discordant})
          </span>
          <span>
            <i className="swatch" style={{ background: "var(--null)", height: "9px", width: "9px", borderRadius: "50%" }} />{" "}
            unresolved in one or both ({q.unresolved})
          </span>
        </div>
      </div>

      <div className="panel">
        <h2 className="panel-title">What the picture says</h2>
        <p className="note" style={{ maxWidth: "70ch" }}>
          Deming slope <span className="mono">{fmtEffect(d.slope, 2).replace("+", "")}</span>{" "}
          <span className="interval">
            [{d.slope_ci_low.toFixed(2)}, {d.slope_ci_high.toFixed(2)}]
          </span>
          . Only <strong>{resolved.length}</strong> of {data.n_pairs.toLocaleString()} ortholog
          groups have an interval that excludes zero in <em>both</em> species at once — and among
          those, <strong>{q.discordant}</strong> point in opposite directions. The honest reading is
          less &ldquo;human and mouse agree&rdquo; than &ldquo;for almost every gene, at least one of
          the two species cannot tell&rdquo;. That is a statement about statistical power across two
          independently pooled corpora, not about biology.
        </p>
      </div>

      <div className="panel">
        <h2 className="panel-title">Genes resolved in both species</h2>
        <div className="table-scroll">
          <table>
            <thead>
              <tr>
                <th>Ortholog</th>
                <th className="num">human g</th>
                <th className="num">mouse g</th>
                <th>Agreement</th>
              </tr>
            </thead>
            <tbody>
              {resolved
                .sort((a, b) => Math.abs(b.g_human) - Math.abs(a.g_human))
                .map((p) => {
                  const agree = p.verdict_human === p.verdict_mouse;
                  return (
                    <tr key={p.ortholog_group}>
                      <td>
                        <button
                          type="button"
                          className="link-button"
                          onClick={() => onPick?.(p.symbol_human)}
                        >
                          {p.symbol_human} / {p.symbol_mouse}
                        </button>
                      </td>
                      <td className={`num ${p.g_human > 0 ? "up" : "down"}`}>{fmtEffect(p.g_human)}</td>
                      <td className={`num ${p.g_mouse > 0 ? "up" : "down"}`}>{fmtEffect(p.g_mouse)}</td>
                      <td className={agree ? "" : "null"}>
                        {agree
                          ? p.verdict_human === "increases"
                            ? "both rise"
                            : "both fall"
                          : "discordant"}
                      </td>
                    </tr>
                  );
                })}
            </tbody>
          </table>
        </div>
        <p className="note" style={{ marginTop: "0.7rem" }}>
          The other {(data.n_pairs - resolved.length).toLocaleString()} ortholog groups are the grey
          cloud: measured in both species, but with an interval crossing zero in at least one, so no
          direction can be claimed.
        </p>
      </div>
    </div>
  );
}
