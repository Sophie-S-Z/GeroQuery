import { useEffect, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import { json } from "../lib/db.js";
import { css } from "../lib/theme.js";
import { fmtInt } from "../lib/format.js";

/**
 * The orientation view — what this corpus is, and the one honest headline the
 * tool had computed and never shown.
 *
 * Of the ~42k genes pooled from the panel, far fewer clear the *prediction*
 * interval (where the next study would land) than clear the confidence interval
 * (where the mean effect is). That gap is the whole thesis of the tool made into
 * a number, and §4.7 flagged that it was being left off the page. It leads here.
 *
 * Everything is read from the shipped artifacts; nothing is computed in the
 * browser beyond laying three precomputed counts end to end.
 */
export default function Overview({ onNavigate }) {
  const host = useRef(null);
  const [meta, setMeta] = useState(null);
  const [crosslayer, setCrosslayer] = useState(null);
  const [crossspecies, setCrossspecies] = useState(null);

  useEffect(() => {
    let live = true;
    json("meta.json").then((d) => live && setMeta(d)).catch(() => {});
    json("crosslayer.json").then((d) => live && setCrosslayer(d)).catch(() => {});
    json("crossspecies.json").then((d) => live && setCrossspecies(d)).catch(() => {});
    return () => {
      live = false;
    };
  }, []);

  const pi = meta?.prediction_interval;

  useEffect(() => {
    if (!host.current || !pi) return undefined;
    const c = css();
    const rows = [
      { stage: "pooled from the panel", n: pi.n_with_pi, kind: "base" },
      { stage: "mean effect excludes zero", n: pi.n_ci_excludes_zero, kind: "ci" },
      { stage: "next study would too", n: pi.n_pi_excludes_zero, kind: "pi" },
    ];
    const chart = Plot.plot({
      width: 680,
      height: 150,
      marginLeft: 196,
      marginRight: 70,
      marginTop: 6,
      marginBottom: 28,
      style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "11px" },
      x: { label: null, grid: true, nice: true },
      y: { label: null, domain: rows.map((r) => r.stage) },
      marks: [
        Plot.barX(rows, {
          y: "stage",
          x: "n",
          fill: (d) => (d.kind === "pi" ? c.up : d.kind === "ci" ? c.inkSoft : c.rule),
          inset: 3,
        }),
        Plot.text(rows, {
          y: "stage",
          x: "n",
          text: (d) => fmtInt(d.n),
          textAnchor: "start",
          dx: 6,
          fill: c.ink,
          fontSize: 11,
          fontWeight: 600,
        }),
      ],
    });
    host.current.replaceChildren(chart);
    return () => chart.remove();
  }, [pi]);

  const findings = [
    crossspecies?.available && {
      view: "crossspecies",
      title: "Human and mouse rarely resolve the same gene",
      body: `Of ${fmtInt(crossspecies.n_pairs)} orthologs measured in both species, only ${
        crossspecies.quadrants.both_up + crossspecies.quadrants.both_down
      } have a direction in both at once — and none disagree.`,
      cue: "See the cross-species scatter →",
    },
    crosslayer?.available && {
      view: "mortality",
      title: "A health state beats nine of ten aging clocks",
      body: `In ${fmtInt(crosslayer.n_subjects)} adults with ${fmtInt(
        crosslayer.n_deaths,
      )} deaths, six routine blood markers predict mortality better than most DNA-methylation clocks.`,
      cue: "See the mortality view →",
    },
    {
      view: "panel",
      title: "One contrast can bend the whole corpus",
      body: "Because the panel is near-balanced, a single study's leave-one-out shift, summed over every gene, shows which contrast moves the published answer most.",
      cue: "See corpus influence →",
    },
  ].filter(Boolean);

  if (!meta) {
    return (
      <div className="panel">
        <div className="skeleton" style={{ width: "60%", height: "2rem" }} />
        <div className="skeleton" style={{ width: "45%" }} />
      </div>
    );
  }

  return (
    <div className="overview">
      <section className="overview-hero">
        <p className="overview-kicker">Measured aging evidence, with intervals</p>
        <h1 className="overview-head">
          {pi ? (
            <>
              <b>{fmtInt(pi.n_pi_excludes_zero)}</b> genes change with age so consistently that the{" "}
              <em>next</em> study would find it too.
            </>
          ) : (
            "Measured aging evidence, with intervals."
          )}
        </h1>
        {pi ? (
          <p className="overview-lede">
            That is the honest headline. {fmtInt(pi.n_ci_excludes_zero)} genes clear the usual
            confidence interval — where the <em>mean</em> effect lies — but only{" "}
            {fmtInt(pi.n_pi_excludes_zero)} clear the prediction interval, where the effect of one
            more study would land. GeroQuery pools {fmtInt(meta.n_genes)} genes from{" "}
            {meta.n_studies} contrasts and shows you which, and why — every number computed in
            Python and read, never re-derived, by this page.
          </p>
        ) : null}
        {pi ? <div className="chart" ref={host} role="img" aria-label="Certainty funnel" /> : null}
      </section>

      <section className="stat-grid">
        {[
          [fmtInt(meta.n_genes), "genes"],
          [fmtInt(meta.n_contrasts_rows), "effect sizes"],
          [`${meta.n_studies} / ${meta.n_series}`, "contrasts / series"],
          [meta.n_tissues, "tissues"],
          [meta.species.length, "species"],
          [crosslayer?.available ? fmtInt(crosslayer.n_deaths) : "—", "linked deaths"],
        ].map(([value, label]) => (
          <div className="stat" key={label}>
            <span className="stat-value">{value}</span>
            <span className="stat-label">{label}</span>
          </div>
        ))}
      </section>

      <section>
        <h2 className="panel-title">Three things a curated database cannot tell you</h2>
        <div className="finding-grid">
          {findings.map((f) => (
            <button
              type="button"
              className="finding-card"
              key={f.view}
              onClick={() => onNavigate?.(f.view)}
            >
              <span className="finding-title">{f.title}</span>
              <span className="finding-body">{f.body}</span>
              <span className="finding-cue">{f.cue}</span>
            </button>
          ))}
        </div>
      </section>

      <section className="panel">
        <h2 className="panel-title">Running it</h2>
        <p className="note" style={{ maxWidth: "74ch" }}>
          This page is static — a React bundle that reads the pooled Parquet in your browser, no
          backend. To rebuild the data from the pinned manifest:{" "}
          <span className="mono">make data &amp;&amp; make idmap</span>, then{" "}
          <span className="mono">python -m geroquery.etl.build_frontend_data</span>. The API and MCP
          server (<span className="mono">uvicorn geroquery.api.app:app</span>,{" "}
          <span className="mono">geroquery-mcp</span>) serve the same estimators. Full instructions
          are in <span className="mono">docs/OVERVIEW.md</span>.
        </p>
        <p className="note" style={{ marginTop: "0.6rem" }}>
          Built {meta.built_on} · manifest {meta.manifest_version} · data version{" "}
          <span className="mono">{meta.data_version}</span>.
        </p>
      </section>
    </div>
  );
}
