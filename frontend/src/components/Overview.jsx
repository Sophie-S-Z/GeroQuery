import { useEffect, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import { json } from "../lib/db.js";
import { css } from "../lib/theme.js";
import { fmtInt, fmtPct } from "../lib/format.js";

/**
 * The orientation view — what this corpus is, and the one honest headline the
 * tool had computed and never shown.
 *
 * Of the ~42k genes pooled from the panel, far fewer clear the prediction
 * interval (where the next study would land) than clear the confidence interval
 * (where the mean effect is). That gap is the whole thesis of the tool made into
 * a number, and §4.7 flagged that it was being left off the page. It leads here.
 *
 * Every number is read from the shipped artifacts; nothing is computed in the
 * browser beyond laying three precomputed counts end to end. The headline count
 * includes genes that rise *and* genes that fall with age, so it is drawn in ink,
 * never in a polarity colour.
 */
export default function Overview({ onNavigate, onPick }) {
  const host = useRef(null);
  const [meta, setMeta] = useState(null);
  const [crosslayer, setCrosslayer] = useState(null);
  const [crossspecies, setCrossspecies] = useState(null);
  const [influence, setInfluence] = useState(null);
  const [studies, setStudies] = useState(null);
  const [curated, setCurated] = useState(null);

  useEffect(() => {
    let live = true;
    const load = (name, set) =>
      json(name)
        .then((d) => live && set(d))
        .catch(() => {});
    load("meta.json", setMeta);
    load("crosslayer.json", setCrosslayer);
    load("crossspecies.json", setCrossspecies);
    load("influence_corpus.json", setInfluence);
    load("studies.json", setStudies);
    load("curated.json", setCurated);
    return () => {
      live = false;
    };
  }, []);

  const pi = meta?.prediction_interval;

  useEffect(() => {
    if (!host.current || !pi) return undefined;
    const c = css();
    // Sequential emphasis, lightest to darkest: the narrowing is the story.
    const rows = [
      { stage: "pooled from the panel", n: pi.n_with_pi, fill: c.rule },
      { stage: "mean effect excludes zero", n: pi.n_ci_excludes_zero, fill: c.inkFaint },
      { stage: "next study would too", n: pi.n_pi_excludes_zero, fill: c.ink },
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
        Plot.barX(rows, { y: "stage", x: "n", fill: "fill", inset: 3 }),
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

  if (!meta) {
    return (
      <div className="panel">
        <div className="skeleton" style={{ width: "60%", height: "2rem" }} />
        <div className="skeleton" style={{ width: "45%" }} />
      </div>
    );
  }

  const top = influence?.[0];
  const topStudy = top && studies?.find((s) => s.study_id === top.study_id);

  const findings = [
    curated?.available && {
      view: "list",
      patch: { genes: curated.genes.map((g) => g.symbol).join(","), species: curated.species },
      title: `${fmtPct(1 - curated.rate_resolved)} of curated aging genes show no direction here`,
      body: `Of ${fmtInt(curated.n_measured)} ${curated.database} genes measured in the panel, ${fmtInt(
        curated.verdicts.increases + curated.verdicts.decreases,
      )} have a pooled direction with age — ${fmtPct(curated.rate_resolved)}, against ${fmtPct(
        curated.background.rate_resolved,
      )} for every ${curated.species} gene. The curated label raises the odds; it does not settle the question.`,
    },
    crossspecies?.available && {
      view: "crossspecies",
      title: "Human and mouse rarely resolve the same gene",
      body: `Of ${fmtInt(crossspecies.n_pairs)} orthologs measured in both species, ${
        crossspecies.quadrants.both_up + crossspecies.quadrants.both_down
      } have a direction in both at once, and ${
        crossspecies.quadrants.discordant
      } of those disagree.`,
    },
    crosslayer?.available && {
      view: "mortality",
      title: "Six blood markers rival the aging clocks",
      body: `In ${fmtInt(crosslayer.n_subjects)} adults followed to ${fmtInt(
        crosslayer.n_deaths,
      )} deaths, a routine blood-chemistry health state predicts mortality better than ${
        crosslayer.clocks.filter((c) => c.c_clock < c.c_dysregulation).length
      } of ${crosslayer.clocks.length} DNA-methylation clocks.`,
    },
    top &&
      topStudy && {
      view: "panel",
      title: "One contrast pulls hardest on the evidence",
      body: `Dropping ${top.study_id.replace(/^GEO:/, "")} (${topStudy.species} ${
        topStudy.tissue
      }) moves the pooled estimate of the ${fmtInt(top.n_genes)} genes it measures by ${top.mean_abs_shift.toFixed(
        3,
      )} on average — more than any other ${topStudy.species} contrast.`,
    },
  ].filter(Boolean);

  return (
    <div className="overview">
      <section className="overview-hero">
        <h1 className="overview-head">
          {pi
            ? `${fmtInt(pi.n_pi_excludes_zero)} of ${fmtInt(pi.n_with_pi)} genes change with age consistently enough that a new study would see it too.`
            : "Measured aging evidence, with intervals."}
        </h1>
        {pi ? (
          <p className="overview-lede">
            {fmtInt(pi.n_ci_excludes_zero)} genes clear the usual confidence interval, which says
            where the average effect lies. Only {fmtInt(pi.n_pi_excludes_zero)} clear the prediction
            interval, which says where one more study would land. GeroQuery re-derives both from{" "}
            {fmtInt(meta.n_contrasts_rows)} published effect sizes across {meta.n_studies}{" "}
            contrasts, {meta.n_tissues} tissues and {meta.species.length} species, and shows the
            evidence behind every gene.
          </p>
        ) : null}

        <form
          className="overview-search"
          role="search"
          onSubmit={(event) => {
            event.preventDefault();
            const symbol = new FormData(event.currentTarget).get("gene")?.toString().trim();
            if (symbol) onPick?.(symbol);
          }}
        >
          <input
            name="gene"
            data-search
            aria-label="Look up a gene"
            placeholder="Look up a gene, e.g. FOXO3"
            autoComplete="off"
            spellCheck="false"
          />
          <span className="kbd" aria-hidden="true" title="Press / to search from anywhere">
            /
          </span>
          <button type="submit">Search</button>
        </form>

        {pi ? <div className="chart" ref={host} role="img" aria-label="Certainty funnel" /> : null}
      </section>

      {findings.length ? (
        <section>
          <h2 className="panel-title">What the corpus says that a curated database cannot</h2>
          <ol className="findings">
            {findings.map((f, i) => (
              <li key={f.view} style={{ "--i": i }}>
                <button type="button" onClick={() => onNavigate?.(f.view, f.patch)}>
                  <span className="finding-title">{f.title}</span>
                  <span className="finding-body">{f.body}</span>
                  <span className="finding-arrow" aria-hidden="true">
                    →
                  </span>
                </button>
              </li>
            ))}
          </ol>
        </section>
      ) : null}

      <section className="panel">
        <h2 className="panel-title">How the numbers are made</h2>
        <p className="note" style={{ maxWidth: "74ch" }}>
          Per-study effects are standardized to Hedges&rsquo; g and pooled by random-effects
          meta-analysis (DerSimonian-Laird with a Hartung-Knapp-Sidik-Jonkman interval, never
          narrower than the conventional one), with a 95% prediction interval for where the next
          study would land. A corpus-wide adaptive-shrinkage prior gives each gene a local false
          sign rate — the probability its direction is wrong. Cross-species agreement uses a Deming
          (errors-in-both-variables) fit; the mortality result is a survey-weighted Cox model with a
          design-based sandwich variance
          {crosslayer?.survival_curves
            ? ", plus weighted Kaplan-Meier curves and restricted mean survival time"
            : ""}
          . Every one of these is computed in Python by the estimator the API serves
          and shipped in the Parquet/JSON; this page renders them and never re-derives them, so the
          website, the API and the MCP tool cannot disagree about what a gene&rsquo;s answer is.
        </p>
      </section>

      <section className="panel">
        <h2 className="panel-title">Running it</h2>
        <p className="note" style={{ maxWidth: "74ch" }}>
          This page is static — a React bundle that reads the pooled Parquet in your browser, no
          backend. To rebuild its data from the pinned manifest:
        </p>
        <pre className="code">
          {"make data && make idmap\npython -m geroquery.etl.build_frontend_data"}
        </pre>
        <p className="note" style={{ maxWidth: "74ch" }}>
          The API (<span className="mono">uvicorn geroquery.api.app:app</span>) and the MCP server (
          <span className="mono">geroquery-mcp</span>) serve the same estimators. Full instructions
          are in <span className="mono">docs/OVERVIEW.md</span>.
        </p>
      </section>
    </div>
  );
}
