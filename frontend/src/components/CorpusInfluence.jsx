import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import { css } from "../lib/theme.js";
import { shortLabel } from "./ForestPlot.jsx";

/**
 * Corpus influence — which of the panel's 32 contrasts bends the whole answer
 * (§4.5). The view with no equivalent in any comparable tool, because no
 * comparable tool re-derives its own numbers.
 *
 * Because the panel is near-balanced — every gene is measured in very nearly the
 * same contrasts — a contrast's leave-one-out shift, summed across all ~42k
 * pooled genes, is a statement about the *panel itself*: drop this one contrast
 * and the corpus moves by this much, in total, across this many genes. That is
 * provenance made measurable, which is this repo's entire argument.
 *
 * The bar is pure magnitude, so it is a single ink mark, not a colour scale;
 * identity (accession, tissue, platform, n) rides in the labels and the table
 * beside it. Everything shown is computed in Python (`corpus_influence`) and
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

    const rows = ranking.map((r) => {
      const meta = byStudy[r.study_id] ?? {};
      return {
        ...r,
        // Compact, still-unique label (two contrasts can share a GDS accession).
        accession: shortLabel(r.study_id),
        tissue: meta.tissue ?? "—",
        species: meta.species ?? "",
        n: meta.sample_size,
      };
    });

    const chart = Plot.plot({
      marginLeft: 150,
      marginRight: 18,
      marginTop: 8,
      marginBottom: 36,
      width: 760,
      height: rows.length * 19 + 54,
      style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "10px" },
      x: {
        label: "total leave-one-out shift, summed over every pooled gene →",
        labelAnchor: "center",
        labelOffset: 32,
        grid: true,
        nice: true,
      },
      y: { label: null, domain: rows.map((r) => r.accession) },
      marks: [
        Plot.ruleX([0], { stroke: c.rule, strokeWidth: 1 }),
        Plot.barX(rows, {
          y: "accession",
          x: "total_abs_shift",
          fill: c.ink,
          inset: 2,
          title: (d) =>
            `${d.accession}\n${d.tissue}${d.species ? ` · ${d.species}` : ""}` +
            `${d.n ? ` · n = ${d.n}` : ""}\n` +
            `total shift ${d.total_abs_shift.toFixed(1)} across ${d.n_genes.toLocaleString()} genes\n` +
            `mean shift per gene ${d.mean_abs_shift.toFixed(4)}`,
        }),
        // The tissue, direct-labelled at the bar end — identity in text, not a
        // second colour encoding on a magnitude bar.
        Plot.text(rows, {
          y: "accession",
          x: "total_abs_shift",
          text: "tissue",
          textAnchor: "start",
          dx: 5,
          fill: c.inkFaint,
          fontSize: 9,
        }),
      ],
    });

    host.current.replaceChildren(chart);
    return () => chart.remove();
  }, [ranking, byStudy]);

  if (!ranking?.length) return null;

  return (
    <div className="chart" ref={host} role="img" aria-label="Corpus influence of each contrast" />
  );
}
