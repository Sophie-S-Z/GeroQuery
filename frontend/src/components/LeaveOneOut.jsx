import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import { css } from "../lib/theme.js";
import { shortLabel } from "./ForestPlot.jsx";

/**
 * The leave-one-out strip — the forest plot's audit companion (§4.4).
 *
 * The forest plot shows where each contrast sits; it cannot show which contrast
 * is *carrying* the pooled result, because a wide row and an influential row
 * look the same. This strip answers that: one row per contrast, in the same
 * order as the forest beside it, each marking where the pooled estimate lands
 * **with that contrast removed**. The full pooled estimate is the heavy
 * reference rule; a dot far from it is a contrast the answer depends on.
 *
 * Every `loo_g` is computed in Python by the same random-effects estimator the
 * pool uses (`geroquery.harmonize.leave_one_out`), shipped in `contrasts.parquet`.
 * This draws it; it does not recompute a pool in the browser.
 *
 * @param {object} props
 * @param {Array<{study_id:string,effect_size:number,loo_g:number,loo_shift:number,tissue:string}>} props.contrasts
 *   Sorted identically to the forest plot (by effect size), each carrying its
 *   leave-one-out columns.
 * @param {{g:number}} props.pooled  The full pooled estimate (the reference).
 * @param {string} props.symbol
 */
export default function LeaveOneOut({ contrasts, pooled, symbol }) {
  const host = useRef(null);

  const usable = (contrasts ?? []).filter((row) => row?.loo_g != null);

  useEffect(() => {
    if (!host.current || !usable.length || !pooled) return undefined;
    const c = css();

    const rows = usable.map((row) => ({
      ...row,
      label: shortLabel(row.study_id),
      // Signed shift: full pooled minus the estimate without this contrast.
      // Reported in Python; recomputed here only as a display fallback.
      shift: row.loo_shift ?? pooled.g - row.loo_g,
    }));

    const values = [...rows.map((r) => r.loo_g), pooled.g, 0];
    const pad = 0.05;
    const domain = [Math.min(...values) - pad, Math.max(...values) + pad];

    const marks = [
      // Zero, faint — so a reader can see whether dropping a contrast ever
      // moves the estimate across the line the whole tool is about.
      Plot.ruleX([0], { stroke: c.rule, strokeWidth: 1 }),
      // The full pooled estimate: the thing every row is compared against.
      Plot.ruleX([pooled.g], { stroke: c.ink, strokeWidth: 1.75 }),
      // The jump: a segment from the full estimate to the leave-one-out one,
      // so influence reads as a length, not a position to be judged by eye.
      Plot.link(rows, {
        y: "label",
        x1: () => pooled.g,
        x2: "loo_g",
        stroke: c.inkFaint,
        strokeWidth: 1.5,
      }),
      Plot.dot(rows, {
        y: "label",
        x: "loo_g",
        fill: c.paper,
        stroke: c.ink,
        strokeWidth: 1.5,
        r: 3.4,
        title: (d) =>
          `${d.study_id.replace(/^GEO:/, "")}\n` +
          `without this contrast: g = ${d.loo_g.toFixed(3)}\n` +
          `shift ${d.shift >= 0 ? "+" : ""}${d.shift.toFixed(3)} from the pooled ${pooled.g.toFixed(3)}`,
      }),
    ];

    const chart = Plot.plot({
      marginLeft: 132,
      marginRight: 16,
      marginTop: 8,
      marginBottom: 34,
      width: 340,
      height: Math.max(180, rows.length * 21 + 54),
      style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "10px" },
      x: {
        domain,
        label: "pooled g with the row removed",
        labelAnchor: "center",
        labelOffset: 30,
        grid: true,
        nice: true,
      },
      y: { label: null, domain: rows.map((r) => r.label) },
      marks,
    });

    host.current.replaceChildren(chart);
    return () => chart.remove();
  }, [usable, pooled, symbol]);

  if (!usable.length || !pooled) return null;

  const biggest = usable.reduce((a, b) =>
    Math.abs(b.loo_shift ?? 0) > Math.abs(a.loo_shift ?? 0) ? b : a,
  );

  return (
    <div className="loo">
      <h3 className="panel-subtitle">Leave-one-out</h3>
      <div className="chart" ref={host} role="img" aria-label={`Leave-one-out audit of ${symbol}`} />
      <p className="note">
        Each dot is the pooled estimate recomputed with that one contrast dropped; the heavy rule is
        the full estimate. Removing{" "}
        <span className="mono">{shortLabel(biggest.study_id)}</span> moves it the most
        {biggest.loo_shift != null
          ? ` (${biggest.loo_shift >= 0 ? "+" : ""}${biggest.loo_shift.toFixed(3)})`
          : ""}
        . A short segment means no single study is carrying the result.
      </p>
    </div>
  );
}
