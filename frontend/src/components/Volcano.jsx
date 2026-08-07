import { useEffect, useMemo, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import { css } from "../lib/theme.js";
import { volcanoPoints } from "../lib/db.js";
import { fmtEffect, fmtInt, fmtLfsr } from "../lib/format.js";

/**
 * Effect size against certainty, for every pooled estimate at once.
 *
 * **Why this view exists.** The landscape histogram bins effect size alone, and
 * that hides the corpus's most important structural fact: the largest effects
 * are also the least certain. A gene can carry a pooled g of +1.15 with an
 * interval that excludes zero and still have close to even odds of pointing the
 * wrong way, because 21,000 genes were measured and some of them were always
 * going to look extreme. You cannot see that in a histogram of g. You can see it
 * here immediately, because those genes sit far right and near the floor.
 *
 * **Why the y-axis is not a p-value.** Two reasons, and the second is a bug.
 *
 * 1. A p-value answers "is this effect exactly zero", which nobody believes of
 *    any gene. The local false sign rate answers "if I claim a direction, how
 *    often am I wrong" — the question a reader actually has.
 * 2. Selecting features by effect size *and* significance — the two thresholds a
 *    volcano plot invites — does not control the false discovery rate.
 *    Benjamini-Hochberg bounds the error rate over all features, not over a
 *    subset carved out with a second filter, and the feature with the largest
 *    estimated effect is a very likely false positive. So there is deliberately
 *    **no "large and significant" selector here**: that convenience is the bug.
 *    The guide lines are reading aids and select nothing.
 *
 * The lfsr is a posterior under a prior fitted across the whole corpus per
 * species, so it already carries the multiplicity that a per-gene interval
 * cannot. That is why far more genes have an interval excluding zero than have
 * a small lfsr, and the gap between those two counts is stated below the chart
 * rather than left for a reader to discover.
 */

/** lfsr = 0 rounds to an infinite ordinate; clamp to the smallest value the
 *  export's five-decimal rounding can actually represent. */
const LFSR_FLOOR = 1e-5;
const GUIDE = 0.05;

export default function Volcano({ species, onPick, onSpecies }) {
  const host = useRef(null);
  const [rows, setRows] = useState([]);
  const [error, setError] = useState(null);

  useEffect(() => {
    let live = true;
    setError(null);
    volcanoPoints(species)
      .then((data) => live && setRows(data))
      .catch((e) => live && setError(e.message));
    return () => {
      live = false;
    };
  }, [species]);

  const points = useMemo(
    () =>
      rows.map((row) => ({
        ...row,
        certainty: -Math.log10(Math.max(row.lfsr, LFSR_FLOOR)),
      })),
    [rows],
  );

  /** The counts the chart is making an argument about, computed once. */
  const summary = useMemo(() => {
    const excludesZero = rows.filter((r) => r.verdict !== "no_evidence").length;
    const confident = rows.filter((r) => r.lfsr <= GUIDE).length;
    const loud = rows.filter((r) => r.verdict !== "no_evidence" && r.lfsr > 0.2).length;
    return { total: rows.length, excludesZero, confident, loud };
  }, [rows]);

  /** Genes worth naming: most certain, not most extreme. Labelling the extremes
   *  would be the winner's curse drawn in text. */
  const labelled = useMemo(
    () =>
      [...points]
        .filter((r) => r.lfsr <= 0.01)
        .sort((a, b) => Math.abs(b.g_shrunk ?? b.g) - Math.abs(a.g_shrunk ?? a.g))
        .slice(0, 12),
    [points],
  );

  useEffect(() => {
    if (!host.current || !points.length) return undefined;
    const c = css();

    const chart = Plot.plot({
      width: 860,
      height: 420,
      marginLeft: 58,
      marginBottom: 42,
      marginRight: 24,
      style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "10px" },
      x: {
        label: "← lower with age    pooled Hedges' g    higher with age →",
        labelAnchor: "center",
        grid: true,
        nice: true,
      },
      y: {
        label: "↑ certainty of direction  −log₁₀(lfsr)",
        grid: true,
        nice: true,
      },
      color: {
        domain: ["decreases", "no_evidence", "increases"],
        range: [c.down, c.null_, c.up],
      },
      marks: [
        // Zero, heavier than the grid, as everywhere else in this app.
        Plot.ruleX([0], { stroke: c.ink, strokeWidth: 1.75 }),
        // The 5% guide. Dashed and unlabelled as a threshold, because it is not
        // one — nothing filters on it.
        Plot.ruleY([-Math.log10(GUIDE)], {
          stroke: c.inkFaint,
          strokeWidth: 1,
          strokeDasharray: "4,4",
        }),
        Plot.dot(points, {
          x: "g",
          y: "certainty",
          fill: "verdict",
          r: 1.5,
          fillOpacity: 0.38,
          title: (d) =>
            `${d.symbol}\ng = ${d.g.toFixed(3)} [${d.ci_low.toFixed(2)}, ${d.ci_high.toFixed(2)}]\n` +
            `shrunk to ${(d.g_shrunk ?? d.g).toFixed(3)}\n` +
            `lfsr = ${d.lfsr.toFixed(4)} — chance the direction is wrong\n` +
            `k = ${d.k} contrasts`,
        }),
        // Named points get a ring so they read above the cloud without a
        // second hue: identity here is position plus a label, not colour.
        Plot.dot(labelled, {
          x: "g",
          y: "certainty",
          fill: "verdict",
          stroke: c.paper,
          strokeWidth: 1.2,
          r: 3.4,
        }),
        Plot.text(labelled, {
          x: "g",
          y: "certainty",
          text: "symbol",
          dy: -9,
          fill: c.ink,
          fontSize: 10,
          fontFamily: c.mono,
        }),
      ],
    });

    host.current.replaceChildren(chart);
    return () => chart.remove();
  }, [points, labelled]);

  return (
    <>
      <div className="panel">
        <h2 className="panel-title">Effect against certainty — {species}</h2>
        {error ? <p className="error">{error}</p> : null}
        <div className="species-switch">
          {["human", "mouse"].map((s) => (
            <button key={s} type="button" aria-pressed={species === s} onClick={() => onSpecies(s)}>
              {s}
            </button>
          ))}
        </div>

        <div
          className="chart"
          ref={host}
          role="img"
          aria-label={`Pooled effect size against local false sign rate for every ${species} gene`}
        />

        <div className="legend">
          <span>
            <i className="swatch" style={{ background: "var(--up)" }} /> interval above zero
          </span>
          <span>
            <i className="swatch" style={{ background: "var(--down)" }} /> interval below zero
          </span>
          <span>
            <i className="swatch" style={{ background: "var(--null)" }} /> interval crosses zero
          </span>
          <span>
            <i className="swatch" style={{ background: "var(--ink-faint)", height: "1px" }} />{" "}
            dashed: lfsr = 0.05
          </span>
        </div>

        {/* `certainty-summary`, not just `note`: this paragraph renders only
            once the data has loaded, so a test selecting the first `.note` on
            the page silently gets the table caption below instead — which is
            what happened, and it passed two of its three assertions against the
            wrong element before failing on the third. */}
        {rows.length ? (
          <p className="note certainty-summary" style={{ marginTop: "0.8rem" }}>
            <strong className="mono">{fmtInt(summary.excludesZero)}</strong> of{" "}
            <strong className="mono">{fmtInt(summary.total)}</strong> genes have an interval that
            excludes zero, but only <strong className="mono">{fmtInt(summary.confident)}</strong>{" "}
            have a local false sign rate at or below 0.05. The two disagree because an interval is a
            statement about one gene and this panel measured{" "}
            <span className="mono">{fmtInt(summary.total)}</span> of them — the lfsr is a posterior
            under a prior fitted across the whole corpus, so it already carries that multiplicity.{" "}
            <strong className="mono">{fmtInt(summary.loud)}</strong> genes have an interval
            excluding zero <em>and</em> more than a one-in-five chance of pointing the wrong way.
            Those are the points sitting far from zero and low on this chart, and they are the ones
            a table ranked by effect size puts at the top.
          </p>
        ) : null}
      </div>

      <div className="panel">
        <h2 className="panel-title">Most certain direction — {species}</h2>
        <div className="table-scroll">
          <table>
            <thead>
              <tr>
                <th>Gene</th>
                <th className="num">g</th>
                <th className="num">shrunk</th>
                <th className="num">lfsr</th>
                <th className="num">k</th>
              </tr>
            </thead>
            <tbody>
              {labelled.map((row) => (
                <tr key={`${row.gene_id}-${row.species}`}>
                  <td>
                    <button
                      type="button"
                      className="link-button mono"
                      onClick={() => onPick(row.symbol)}
                    >
                      {row.symbol}
                    </button>
                  </td>
                  <td className={`num ${row.g > 0 ? "up" : "down"}`}>{fmtEffect(row.g)}</td>
                  <td className="num">{fmtEffect(row.g_shrunk)}</td>
                  <td className="num">{fmtLfsr(row.lfsr)}</td>
                  <td className="num">{row.k}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <p className="note" style={{ marginTop: "0.7rem" }}>
          Ranked by shrunken effect among genes with lfsr ≤ 0.01 — large <em>and</em> well measured.
          The <span className="mono">shrunk</span> column is the posterior mean once the corpus
          prior is applied; the gap between it and <span className="mono">g</span> is how much of
          the raw estimate was noise. It is largest exactly where the raw estimate is largest.
        </p>
      </div>
    </>
  );
}
