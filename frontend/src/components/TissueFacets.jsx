import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import { css } from "../lib/theme.js";
import { shortLabel } from "./ForestPlot.jsx";

/**
 * Tissue small-multiples for one gene (§4.6).
 *
 * The pooled number hides that, for a human gene, most of the evidence is often
 * one tissue — CDKN2A's 14 human contrasts are 10 skeletal muscle and one each
 * of three others. A single pooled diamond cannot show that; small multiples
 * can. One mini-forest per tissue, a **shared x scale** so magnitudes compare
 * across panels, sorted by contrast count so the load-bearing tissue is first.
 *
 * Small multiples, not a colour-per-tissue scale: mouse spans nine tissues,
 * above the eight-hue ceiling, and a ninth categorical hue is never generated
 * (§7). Facetting is both the rule and the clearer answer.
 *
 * The lean panels are drawn, not dropped: a tissue with one contrast *is* the
 * limitation section made visible, and hiding it would hide the point. Nothing
 * is pooled per tissue here — these are the raw per-study effects already in
 * `contrasts.parquet`, with the gene's overall pooled estimate as a reference.
 */
export default function TissueFacets({ contrasts, pooled, symbol }) {
  const host = useRef(null);

  useEffect(() => {
    if (!host.current || !contrasts?.length) return undefined;
    const c = css();
    const colour = { increases: c.up, decreases: c.down, no_evidence: c.null_ };

    const rows = contrasts.map((row) => {
      const half = 1.96 * (row.standard_error ?? 0);
      const low = row.effect_size - half;
      const high = row.effect_size + half;
      return {
        ...row,
        label: shortLabel(row.study_id),
        low,
        high,
        tissue: row.tissue ?? "unspecified",
        dir: low > 0 ? "increases" : high < 0 ? "decreases" : "no_evidence",
      };
    });

    // One shared x domain across every tissue, so a dot's position means the
    // same thing in every panel.
    const allX = [...rows.flatMap((r) => [r.low, r.high]), 0, ...(pooled ? [pooled.g] : [])];
    const pad = 0.15;
    const domain = [Math.min(...allX) - pad, Math.max(...allX) + pad];

    // Tissues ordered by contrast count, descending — the load-bearing one first.
    const counts = new Map();
    for (const r of rows) counts.set(r.tissue, (counts.get(r.tissue) ?? 0) + 1);
    const tissues = [...counts.keys()].sort((a, b) => counts.get(b) - counts.get(a) || a.localeCompare(b));

    const charts = tissues.map((tissue) => {
      const here = rows.filter((r) => r.tissue === tissue).sort((a, b) => a.effect_size - b.effect_size);
      const marks = [
        Plot.ruleX([0], { stroke: c.ink, strokeWidth: 1.75 }),
        ...(pooled
          ? [Plot.ruleX([pooled.g], { stroke: c.inkFaint, strokeWidth: 1, strokeDasharray: "3,3" })]
          : []),
        Plot.ruleY(here, {
          y: "label",
          x1: "low",
          x2: "high",
          stroke: (d) => colour[d.dir],
          strokeWidth: 1.4,
          strokeLinecap: "round",
        }),
        Plot.dot(here, {
          y: "label",
          x: "effect_size",
          fill: (d) => (d.dir === "no_evidence" ? c.paper : colour[d.dir]),
          stroke: (d) => colour[d.dir],
          strokeWidth: 1.4,
          r: 3.2,
          title: (d) =>
            `${d.study_id.replace(/^GEO:/, "")}\n${tissue}\n` +
            `g = ${d.effect_size.toFixed(3)} [${d.low.toFixed(2)}, ${d.high.toFixed(2)}]`,
        }),
      ];
      const chart = Plot.plot({
        marginLeft: 118,
        marginRight: 14,
        marginTop: 26,
        marginBottom: 30,
        width: 440,
        height: here.length * 20 + 58,
        style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "10px" },
        x: { domain, label: null, grid: true, nice: true },
        y: { label: null, domain: here.map((r) => r.label) },
        marks,
      });
      const panel = document.createElement("figure");
      panel.className = "facet";
      const caption = document.createElement("figcaption");
      caption.className = "facet-title";
      caption.textContent = `${tissue} · ${counts.get(tissue)} contrast${counts.get(tissue) > 1 ? "s" : ""}`;
      panel.append(caption, chart);
      return panel;
    });

    host.current.replaceChildren(...charts);
    return () => host.current?.replaceChildren();
  }, [contrasts, pooled, symbol]);

  if (!contrasts?.length) return null;

  return (
    <>
      <div className="facet-grid" ref={host} role="img" aria-label={`Per-tissue effects for ${symbol}`} />
      <p className="note">
        Each panel is one tissue; the dashed rule is the gene&rsquo;s overall pooled estimate, the
        heavy rule is zero. A panel with a single contrast is kept, not dropped — the imbalance
        between tissues is itself what a single pooled number hides.
      </p>
    </>
  );
}
