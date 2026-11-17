import { useEffect, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import { json, knockdownsScoredTwice, replicatingKnockdowns } from "../lib/db.js";
import { css } from "../lib/theme.js";

/**
 * The perturbation layer — a separate surface, on purpose.
 *
 * Every other view in this app reports what was measured about age. This one
 * reports what happens when a gene is switched off in a cell line, which is a
 * different experiment answering a different question, so it is kept apart and
 * never folded into a gene's aging verdict.
 *
 * The view leads with what the data does not support. Two cell lines agree on
 * the direction of a knockdown's aging alignment for 54.6% of genes — a coin
 * flip — so the scatter of one screen against the other *is* the finding, and
 * the short list of genes that do replicate sits beneath it as leads rather
 * than results. Every number shown is computed in
 * `geroquery/etl/build_perturbation.py`; nothing is derived here.
 */
export default function Perturbation({ onPick }) {
  const host = useRef(null);
  const [meta, setMeta] = useState(null);
  const [both, setBoth] = useState(null);
  const [hits, setHits] = useState(null);
  const [error, setError] = useState(null);

  useEffect(() => {
    let live = true;
    Promise.all([json("meta.json"), knockdownsScoredTwice(), replicatingKnockdowns()])
      .then(([m, scored, replicating]) => {
        if (!live) return;
        setMeta(m.perturbation ?? null);
        setBoth(scored);
        setHits(replicating);
      })
      .catch((e) => live && setError(e.message));
    return () => {
      live = false;
    };
  }, []);

  useEffect(() => {
    if (!host.current || !both?.length) return undefined;
    const c = css();
    const chart = Plot.plot({
      width: 620,
      height: 620,
      marginLeft: 56,
      marginBottom: 48,
      marginTop: 30,
      marginRight: 16,
      style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "11px" },
      x: { label: "RPE1 adjusted score →", grid: true, nice: true },
      y: { label: "↑ K562 adjusted score", grid: true, nice: true },
      marks: [
        Plot.ruleX([0], { stroke: c.rule, strokeWidth: 1 }),
        Plot.ruleY([0], { stroke: c.rule, strokeWidth: 1 }),
        Plot.dot(both, {
          x: "rpe1_z",
          y: "k562_z",
          r: 2,
          fill: (d) => (d.replicates ? (d.rpe1_z > 0 ? c.up : c.down) : c.null_),
          fillOpacity: (d) => (d.replicates ? 0.95 : 0.3),
          title: (d) =>
            `${d.symbol}\nRPE1 ${d.rpe1_z.toFixed(2)}  K562 ${d.k562_z.toFixed(2)}` +
            (d.replicates ? "\nextreme and agreeing in both screens" : ""),
        }),
      ],
    });
    host.current.replaceChildren(chart);
    return () => chart.remove();
  }, [both]);

  if (error) return <p className="note">Could not load the perturbation layer: {error}</p>;
  if (!meta || !both || !hits) return <p className="note">Loading…</p>;

  const [rpe1, k562] = meta.screens;
  const coverage = meta.coverage;
  const untested = coverage.resolved_absent + coverage.resolved_readout_only;

  return (
    <div className="stack">
      <div className="panel">
        <h2 className="panel-title">What happens when the gene is switched off</h2>
        <p className="note" style={{ maxWidth: "74ch" }}>
          The panel elsewhere in this site is observational: it reads how expression differs between
          young and old tissue, which says nothing about what a gene <em>does</em>. This layer adds
          CRISPRi knockdown screens, where one gene at a time was switched off in a human cell line
          and the whole transcriptome read out, and asks whether the genes that move afterwards move
          the way they move with age.
        </p>
        <p className="note" style={{ maxWidth: "74ch" }}>
          <strong>Neither cell line is aged tissue.</strong> RPE1 is retinal pigment epithelium and
          K562 is a leukaemia line; both are immortalized dividing cells. Nothing here measures
          aging, and no result on this page changes any gene&rsquo;s verdict. These are leads.
        </p>
      </div>

      <div className="panel">
        <h2 className="panel-title">Which aging genes can be asked at all</h2>
        <div className="table-scroll">
          <table>
            <thead>
              <tr>
                <th>Of the {coverage.resolved_genes.toLocaleString()} genes with a direction</th>
                <th className="num">genes</th>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td>Switched off, and the cell measurably responded</td>
                <td className="num">{coverage.resolved_scored.toLocaleString()}</td>
              </tr>
              <tr>
                <td>Switched off, and nothing measurable happened</td>
                <td className="num">{coverage.resolved_no_effect.toLocaleString()}</td>
              </tr>
              <tr>
                <td>Measured as a readout gene, never itself switched off</td>
                <td className="num">{coverage.resolved_readout_only.toLocaleString()}</td>
              </tr>
              <tr>
                <td>Absent from both screens</td>
                <td className="num">{coverage.resolved_absent.toLocaleString()}</td>
              </tr>
            </tbody>
          </table>
        </div>
        <p className="note" style={{ marginTop: "0.7rem", maxWidth: "74ch" }}>
          For <strong>{untested.toLocaleString()}</strong> of them no one has published a knockdown
          of this kind, so this layer has nothing to say about them — which is the most durable
          thing on this page.
        </p>
      </div>

      <div className="panel">
        <h2 className="panel-title">The two cell lines mostly disagree</h2>
        <div ref={host} className="chart" />
        <p className="note" style={{ maxWidth: "74ch" }}>
          Each point is one of the {meta.n_scored_in_both.toLocaleString()} genes switched off in
          both screens, scored against the pooled aging direction. The score is adjusted for how
          hard the knockdown hit the cell, because the raw version correlates{" "}
          <span className="mono">+{rpe1.magnitude_confound.toFixed(2)}</span> and{" "}
          <span className="mono">+{k562.magnitude_confound.toFixed(2)}</span> with the size of the
          response: damaging a cell makes it look older whatever you damaged, and ranking on that
          would rank damage and call it biology.
        </p>
        <p className="note" style={{ maxWidth: "74ch" }}>
          The two screens agree on the sign for{" "}
          <strong>{(meta.cross_screen_sign_agreement * 100).toFixed(1)}%</strong> of genes
          (ρ = <span className="mono">{meta.cross_screen_spearman.toFixed(2)}</span>), and no gene
          reaches a significant alignment in either screen after correction. Only the extremes
          reproduce, so a gene is marked below only when <em>both</em> screens place it in their
          most extreme {Math.round(meta.extreme_quantile * 100)}% and agree on the direction.
        </p>
      </div>

      <div className="panel">
        <h2 className="panel-title">
          Genes whose knockdown replicates across both cell lines ({meta.n_replicating})
        </h2>
        <div className="table-scroll">
          <table>
            <thead>
              <tr>
                <th>Gene</th>
                <th className="num">RPE1</th>
                <th className="num">K562</th>
                <th>Direction</th>
              </tr>
            </thead>
            <tbody>
              {hits.map((row) => (
                <tr key={row.gene_id}>
                  <td>
                    <button
                      type="button"
                      className="link-button"
                      onClick={() => onPick?.(row.symbol)}
                    >
                      {row.symbol}
                    </button>
                  </td>
                  <td className="num">{row.rpe1_z.toFixed(2)}</td>
                  <td className="num">{row.k562_z.toFixed(2)}</td>
                  <td className={row.direction === "toward_aging" ? "up" : "down"}>
                    {row.direction === "toward_aging"
                      ? "switching it off looks more like age"
                      : "switching it off looks less like age"}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <p className="note" style={{ marginTop: "0.7rem", maxWidth: "74ch" }}>
          {meta.n_replicating} of {meta.n_scored_in_both.toLocaleString()} genes scored twice. Read
          these as hypotheses worth a real experiment, not as evidence that switching the gene off
          would do anything to an aging animal.
        </p>
      </div>

      <div className="panel">
        <h2 className="panel-title">Provenance</h2>
        <p className="note" style={{ maxWidth: "74ch" }}>
          {meta.release}, {meta.license}. {rpe1.n_perturbations.toLocaleString()} knockdowns in{" "}
          {rpe1.cell_line} and {k562.n_perturbations.toLocaleString()} in {k562.cell_line}, of which{" "}
          {(rpe1.n_scored + k562.n_scored).toLocaleString()} passed the screens&rsquo; own test that
          anything changed. Each knockdown lowers its own target in{" "}
          {(rpe1.on_target_fraction_down * 100).toFixed(1)}% and{" "}
          {(k562.on_target_fraction_down * 100).toFixed(1)}% of cases respectively, which is the
          control that makes the negative results above readable as measurements. Method and the
          full working: <span className="mono">docs/PERTURBATION.md</span>.
        </p>
        <p className="note" style={{ maxWidth: "74ch" }}>{meta.attribution}</p>
      </div>
    </div>
  );
}
