import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import { css } from "../lib/theme.js";
import { withShortLabels } from "./ForestPlot.jsx";

// Bars per species. Past a few dozen a bar chart stops being read and starts
// being scrolled; the table below carries the rest.
const TOP = 20;

/**
 * Corpus influence — which of the panel's contrasts bends the whole answer
 * (§4.5). The view with no equivalent in any comparable tool, because no
 * comparable tool re-derives its own numbers.
 *
 * Each bar is a contrast's mean leave-one-out shift over the genes it measures:
 * drop this one contrast and the pooled estimates that depend on it move this
 * far on average. The mean, not the total, because coverage varies twentyfold
 * between arrays and a total mostly ranks array size. That is provenance made
 * measurable, which is this repo's entire argument.
 *
 * The bar is pure magnitude, so it is a single ink mark, not a colour scale;
 * identity (accession, tissue, platform, n) rides in the labels and the table
 * beside it. Everything shown is computed in Python (`add_influence`) and
 * shipped in `influence_corpus.json`.
 *
 * @param {object} props
 * @param {Array<{study_id:string,total_abs_shift:number,mean_abs_shift:number,n_genes:number}>} props.ranking
 * @param {Record<string, {tissue?:string, species?:string, sample_size?:number}>} props.byStudy
 */
export default function CorpusInfluence({ ranking, byStudy }) {
  const host = useRef(null);

  useEffect(() => {
    if (!host.current || !ranking?.length) return undefined;
    const c = css();

    const rows = withShortLabels(ranking).map((r) => {
      const meta = byStudy[r.study_id] ?? {};
      return {
        ...r,
        accession: r.label,
        tissue: meta.tissue ?? "—",
        species: meta.species ?? "",
        n: meta.sample_size,
      };
    });

    // One panel per species on a shared x scale. A mouse gene is pooled from
    // fewer, noisier contrasts than a human one, so dropping a mouse contrast
    // moves it further by construction: bars compare only within a species.
    const max = Math.max(...rows.map((r) => r.mean_abs_shift));
    const node = host.current;
    const panels = [...new Set(rows.map((r) => r.species))].sort().map((species) => {
      // The ranking arrives sorted by mean shift; the chart shows the head of
      // it and the provenance table below lists every contrast's shift.
      const all = rows.filter((r) => r.species === species);
      const here = all.slice(0, TOP);
      const chart = Plot.plot({
        marginLeft: 150,
        // Room for the tissue label written past the longest bar.
        marginRight: 96,
        marginTop: 8,
        marginBottom: 36,
        width: 760,
        height: here.length * 19 + 54,
        style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "10px" },
        x: {
          domain: [0, max],
          label: "mean leave-one-out shift per gene measured →",
          labelAnchor: "center",
          labelOffset: 32,
          grid: true,
          nice: true,
        },
        y: { label: null, domain: here.map((r) => r.accession) },
        marks: [
          Plot.ruleX([0], { stroke: c.rule, strokeWidth: 1 }),
          Plot.barX(here, {
            y: "accession",
            x: "mean_abs_shift",
            fill: c.ink,
            inset: 2,
            title: (d) =>
              `${d.accession}\n${d.tissue}${d.n ? ` · n = ${d.n}` : ""}\n` +
              `mean shift per gene ${d.mean_abs_shift.toFixed(4)} across ${d.n_genes.toLocaleString()} genes\n` +
              `total ${d.total_abs_shift.toFixed(1)}`,
          }),
          // The tissue, direct-labelled at the bar end — identity in text, not a
          // second colour encoding on a magnitude bar.
          Plot.text(here, {
            y: "accession",
            x: "mean_abs_shift",
            text: "tissue",
            textAnchor: "start",
            dx: 5,
            fill: c.inkFaint,
            fontSize: 9,
          }),
        ],
      });
      const figure = document.createElement("figure");
      figure.className = "facet";
      const caption = document.createElement("figcaption");
      caption.className = "facet-title";
      caption.textContent =
        here.length < all.length
          ? `${species || "unspecified"} · top ${here.length} of ${all.length} contrasts`
          : `${species || "unspecified"} · ${all.length} contrasts`;
      figure.append(caption, chart);
      return figure;
    });

    node.replaceChildren(...panels);
    return () => node.replaceChildren();
  }, [ranking, byStudy]);

  if (!ranking?.length) return null;

  return (
    <div className="chart" ref={host} role="img" aria-label="Corpus influence of each contrast" />
  );
}
