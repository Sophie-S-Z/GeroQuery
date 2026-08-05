import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import { css } from "../lib/theme.js";

/**
 * The forest plot. One row per contrast, each a 95% interval around a
 * standardized effect, with the pooled estimate as a diamond beneath.
 *
 * The zero rule is drawn heavier than the axis and every other line on the
 * page. That is the design thesis made literal: a curated database can only
 * say yes or no, and the thing this tool adds is an interval that can straddle
 * that line and say "we cannot tell".
 *
 */

/**
 * A short, *unique* axis label for one contrast.
 *
 * Study ids are `GEO:GDS3976` or `GEO:GDS3976:hematopoietic_stem_cell_hsc` —
 * one GDS can yield more than one contrast. Stripping to the bare accession
 * would give those two rows the same label, and Plot's ordinal y-scale would
 * silently collapse them into a single row: two distinct measurements drawn as
 * one. So the qualifier is kept, abbreviated rather than dropped, and the full
 * id stays in the tooltip.
 */
function shortLabel(studyId) {
  const [, accession, qualifier] = studyId.split(":");
  if (!qualifier) return accession;
  const compact = qualifier
    .split("_")
    // Initials for all but the last word: hematopoietic_stem_cell_hsc -> h.s.c.hsc
    .map((word, index, all) => (index === all.length - 1 ? word : word[0]))
    .join(".");
  return `${accession} ${compact.length > 14 ? `${compact.slice(0, 13)}…` : compact}`;
}

/**
 * @param {object} props
 * @param {Array<{study_id:string,effect_size:number,standard_error:number,tissue:string}>} props.contrasts
 * @param {{g:number,ci_low:number,ci_high:number,k:number}|null} props.pooled
 * @param {string} props.symbol
 */
export default function ForestPlot({ contrasts, pooled, symbol }) {
  const host = useRef(null);

  useEffect(() => {
    if (!host.current || !contrasts?.length) return undefined;
    const c = css();

    const rows = contrasts.map((row) => {
      const half = 1.96 * (row.standard_error ?? 0);
      return {
        ...row,
        label: shortLabel(row.study_id),
        low: row.effect_size - half,
        high: row.effect_size + half,
        // Colour by this contrast's own interval, not by the pooled direction.
        // A forest plot whose rows are all tinted by the pooled verdict hides
        // the disagreement that produced it — and for CDKN2A the disagreement
        // is the finding: six up, eight down.
        dir:
          row.effect_size - half > 0
            ? "increases"
            : row.effect_size + half < 0
              ? "decreases"
              : "no_evidence",
      };
    });

    const values = [
      ...rows.flatMap((r) => [r.low, r.high]),
      ...(pooled ? [pooled.ci_low, pooled.ci_high] : []),
      0,
    ];
    const pad = 0.15;
    const domain = [Math.min(...values) - pad, Math.max(...values) + pad];

    const colour = {
      increases: c.up,
      decreases: c.down,
      no_evidence: c.null_,
    };

    const marks = [
      // The zero rule, first so everything else sits above it.
      Plot.ruleX([0], { stroke: c.ink, strokeWidth: 1.75 }),
      Plot.ruleY(rows, {
        y: "label",
        x1: "low",
        x2: "high",
        stroke: (d) => colour[d.dir],
        strokeWidth: 1.5,
        strokeLinecap: "round",
      }),
      Plot.dot(rows, {
        y: "label",
        x: "effect_size",
        fill: (d) => (d.dir === "no_evidence" ? c.paper : colour[d.dir]),
        stroke: (d) => colour[d.dir],
        strokeWidth: 1.5,
        r: 3.6,
        // The full study id, not the abbreviated axis label: the tooltip is
        // where the qualifier that got shortened has to remain readable.
        title: (d) =>
          `${d.study_id.replace(/^GEO:/, "")}\n${d.tissue ?? "tissue unspecified"}\n` +
          `g = ${d.effect_size.toFixed(3)} [${d.low.toFixed(2)}, ${d.high.toFixed(2)}]` +
          (d.age_range ? `\n${d.age_range}` : ""),
      }),
    ];

    if (pooled) {
      const pooledRow = [{ label: "POOLED", ...pooled }];
      marks.push(
        Plot.ruleY(pooledRow, {
          y: "label",
          x1: "ci_low",
          x2: "ci_high",
          stroke: c.ink,
          strokeWidth: 2.5,
        }),
        // A diamond, the forest-plot convention for a summary estimate — so it
        // reads as a different kind of thing from the studies above it.
        Plot.dot(pooledRow, {
          y: "label",
          x: "g",
          symbol: "diamond",
          fill: c.ink,
          r: 6,
          title: (d) =>
            `Pooled (random effects, DerSimonian-Laird)\n` +
            `g = ${d.g} [${d.ci_low}, ${d.ci_high}]\nk = ${d.k} contrasts`,
        }),
      );
    }

    const chart = Plot.plot({
      marginLeft: 132,
      marginRight: 18,
      marginTop: 8,
      marginBottom: 34,
      width: 720,
      height: Math.max(180, rows.length * 21 + (pooled ? 46 : 0) + 54),
      style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "10px" },
      x: {
        domain,
        label: "← lower with age    Hedges' g    higher with age →",
        labelAnchor: "center",
        labelOffset: 30,
        grid: true,
        nice: true,
      },
      y: {
        label: null,
        domain: [...rows.map((r) => r.label), ...(pooled ? ["POOLED"] : [])],
      },
      marks,
    });

    host.current.replaceChildren(chart);
    return () => chart.remove();
  }, [contrasts, pooled, symbol]);

  if (!contrasts?.length) {
    return (
      <p className="note">
        No per-study contrasts for this gene and species.
      </p>
    );
  }

  return (
    <>
      <div className="chart" ref={host} role="img" aria-label={`Forest plot of ${symbol}`} />
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
      </div>
    </>
  );
}
