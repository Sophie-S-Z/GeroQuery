import { useEffect, useMemo, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import { json, panelSummary, pooledForSymbols } from "../lib/db.js";
import { downloadCsv } from "../lib/export.js";
import {
  VERDICT_LABEL,
  fmtEffect,
  fmtInt,
  fmtInterval,
  fmtLfsr,
  fmtPct,
} from "../lib/format.js";
import { css } from "../lib/theme.js";

const MAX_GENES = 2000;

// Spreadsheets silently turn symbols such as SEPT1 and MARCH1 into dates
// (Ziemann et al., Genome Biology 2016). A pasted "1-Sep" is a gene lost before
// it reached this page, so it is named as such, not filed as "not measured".
const MONTHS = "jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec";
const DATE_LIKE = new RegExp(`^(\\d{1,2}-(${MONTHS})|(${MONTHS})-\\d{1,2})$`, "i");

/** Unique symbols in the order given, however the list was separated. */
function parseList(text) {
  const seen = new Set();
  const out = [];
  for (const token of text.split(/[\s,;]+/)) {
    const key = token.toUpperCase();
    if (!token || seen.has(key)) continue;
    seen.add(key);
    out.push(token);
  }
  return out.slice(0, MAX_GENES);
}

const resolvesNextStudy = (row) => row.pi_low > 0 || row.pi_high < 0;

/**
 * Check a gene list against the pooled evidence.
 *
 * The question a researcher actually arrives with is rarely one gene: it is a
 * hit list from their own experiment, a pathway, or a curated set. This looks
 * every symbol up, counts the verdicts the pipeline already computed, and sets
 * the list's rate beside the rate for every gene of the species. Counting is
 * rendering; no estimate is made here, and deliberately no enrichment p-value —
 * a list chosen by any criterion is not a random draw from the panel.
 */
export default function GeneList({ genes, species, onUpdate, onPick }) {
  const symbols = useMemo(() => parseList(genes), [genes]);
  const [rows, setRows] = useState(null);
  const [summary, setSummary] = useState(null);
  const [curated, setCurated] = useState(null);
  const [meta, setMeta] = useState(null);
  const [error, setError] = useState(null);
  const host = useRef(null);

  useEffect(() => {
    let live = true;
    panelSummary()
      .then((s) => live && setSummary(s))
      .catch((e) => live && setError(e.message));
    json("curated.json")
      .then((c) => live && setCurated(c))
      .catch(() => {});
    json("meta.json")
      .then((m) => live && setMeta(m))
      .catch(() => {});
    return () => {
      live = false;
    };
  }, []);

  useEffect(() => {
    let live = true;
    setRows(null);
    if (!symbols.length) return undefined;
    pooledForSymbols(symbols, species)
      .then((r) => live && setRows(r))
      .catch((e) => live && setError(e.message));
    return () => {
      live = false;
    };
  }, [symbols, species]);

  const evidence = useMemo(
    () => new Map((curated?.available ? curated.genes : []).map((g) => [g.symbol, g.evidence])),
    [curated],
  );

  const result = useMemo(() => {
    if (!rows) return null;
    const byKey = new Map();
    for (const row of rows) {
      const key = row.symbol.toUpperCase();
      if (!byKey.has(key)) byKey.set(key, row);
    }
    const found = [...byKey.values()].sort((a, b) => (a.lfsr ?? 1) - (b.lfsr ?? 1));
    const missing = symbols.filter((s) => !byKey.has(s.toUpperCase()));
    const counts = { increases: 0, decreases: 0, no_evidence: 0 };
    for (const row of found) counts[row.verdict] += 1;
    return {
      found,
      missing,
      counts,
      resolved: counts.increases + counts.decreases,
      nextStudy: found.filter(resolvesNextStudy).length,
      dates: missing.filter((s) => DATE_LIKE.test(s)),
    };
  }, [rows, symbols]);

  const background = summary?.find((s) => s.species === species);

  useEffect(() => {
    if (!host.current || !result?.found.length || !background) return undefined;
    const c = css();
    const groups = [
      ["this list", result.found.length, result.counts.increases, result.counts.decreases],
      [`all ${species} genes`, background.n_pooled, background.n_up, background.n_down],
    ];
    const data = groups.flatMap(([group, n, up, down]) => [
      { group, verdict: "increases", share: up / n },
      { group, verdict: "decreases", share: down / n },
      { group, verdict: "no_evidence", share: (n - up - down) / n },
    ]);
    const chart = Plot.plot({
      width: 640,
      height: 92,
      marginLeft: 118,
      marginRight: 14,
      marginTop: 4,
      marginBottom: 24,
      style: { background: "transparent", color: c.inkFaint, fontFamily: c.mono, fontSize: "11px" },
      x: { domain: [0, 1], tickFormat: "%", label: null, grid: true },
      y: { label: null, domain: groups.map(([group]) => group) },
      color: { domain: ["increases", "decreases", "no_evidence"], range: [c.up, c.down, c.null_] },
      marks: [
        Plot.barX(
          data,
          Plot.stackX({
            y: "group",
            x: "share",
            fill: "verdict",
            inset: 2,
            title: (d) => `${d.group} · ${VERDICT_LABEL[d.verdict]}: ${fmtPct(d.share)}`,
          }),
        ),
      ],
    });
    host.current.replaceChildren(chart);
    return () => chart.remove();
  }, [result, background, species]);

  const loadGenAge = () =>
    onUpdate({ genes: curated.genes.map((g) => g.symbol).join(","), species: "human" });

  const download = () =>
    downloadCsv(
      `geroquery-gene-list-${species}-${meta?.data_version ?? "data"}.csv`,
      ["symbol", "species", "verdict", "g", "ci_low", "ci_high", "pi_low", "pi_high", "k", "n_contrasts", "lfsr", "genage_evidence"],
      result.found.map((row) => ({ ...row, genage_evidence: evidence.get(row.symbol) })),
    );

  return (
    <div className="split">
      <form
        className="rail"
        onSubmit={(event) => {
          event.preventDefault();
          const text = new FormData(event.currentTarget).get("genes")?.toString() ?? "";
          onUpdate({ genes: parseList(text).join(",") });
        }}
      >
        <label htmlFor="gene-list" className="panel-title">
          Gene list
        </label>
        <textarea
          id="gene-list"
          key={genes}
          name="genes"
          className="search-field list-field"
          rows={10}
          defaultValue={symbols.join("\n")}
          placeholder={"FOXO3\nCDKN2A\nIGF1R"}
          spellCheck="false"
          autoComplete="off"
        />
        <div className="species-switch" style={{ margin: "0.6rem 0 0" }}>
          {(summary ?? []).map((s) => (
            <button
              key={s.species}
              type="button"
              aria-pressed={s.species === species}
              onClick={() => onUpdate({ species: s.species })}
            >
              {s.species}
            </button>
          ))}
        </div>
        <div className="list-actions">
          <button type="submit" className="button-primary">
            Check list
          </button>
          {curated?.available && symbols.length ? (
            <button type="button" className="button-quiet" onClick={loadGenAge}>
              Load GenAge ({fmtInt(curated.n_listed)})
            </button>
          ) : null}
        </div>
        <p className="rail-hint">
          Separate symbols by lines, spaces or commas — up to {fmtInt(MAX_GENES)}. Matching
          ignores case. The list is checked in your browser and kept in the page address, so the
          link reproduces the check.
        </p>
      </form>

      <div>
        {error ? <p className="error">{error}</p> : null}

        {!symbols.length ? (
          <div className="state">
            <p className="state-title">Check a gene list against the evidence</p>
            <p className="note" style={{ maxWidth: "52ch", margin: "0 auto" }}>
              Paste the genes from your own work — a differential-expression hit list, a pathway,
              a candidate panel — to see which have a pooled direction with age, and how often
              that happens compared with every gene measured here.
            </p>
            {curated?.available ? (
              <p style={{ marginTop: "1rem" }}>
                <button type="button" className="button-quiet" onClick={loadGenAge}>
                  Try GenAge&rsquo;s {fmtInt(curated.n_listed)} curated human aging genes
                </button>
              </p>
            ) : null}
          </div>
        ) : null}

        {symbols.length && !result ? (
          <div className="panel" aria-busy="true">
            <div className="skeleton" style={{ width: "40%", height: "2rem" }} />
            <div className="skeleton" style={{ width: "70%" }} />
          </div>
        ) : null}

        {result && background ? (
          <>
            <h1 className="section-title">
              {fmtInt(result.found.length)} of {fmtInt(symbols.length)} genes measured
            </h1>
            <p className="overview-lede list-lede">
              {fmtInt(result.resolved)} have a pooled direction with age (
              {fmtPct(result.found.length && result.resolved / result.found.length)}), against{" "}
              {fmtPct((background.n_up + background.n_down) / background.n_pooled)} of all{" "}
              {fmtInt(background.n_pooled)} {species} genes. {fmtInt(result.nextStudy)} clear the
              prediction interval, so a new study would be expected to see them too (
              {fmtPct(result.found.length && result.nextStudy / result.found.length)}, against{" "}
              {fmtPct(background.n_pi / background.n_pooled)}).
            </p>

            {result.found.length ? (
              <div className="panel">
                <h2 className="panel-title">This list beside every {species} gene</h2>
                <div
                  className="chart"
                  ref={host}
                  role="img"
                  aria-label={`Share of the list and of all ${species} genes that rise, fall, or show no direction with age`}
                />
                <div className="legend">
                  <span>
                    <i className="swatch" style={{ background: "var(--up)" }} /> rises with age
                  </span>
                  <span>
                    <i className="swatch" style={{ background: "var(--down)" }} /> falls with age
                  </span>
                  <span>
                    <i className="swatch" style={{ background: "var(--null)" }} /> no direction
                  </span>
                </div>
                <p className="note" style={{ marginTop: "0.7rem", maxWidth: "74ch" }}>
                  A lookup and two rates, not an enrichment test. A list chosen by any criterion —
                  significance in your own data, membership of a curated database — is not a
                  random draw from these genes, so no p-value against it would have a valid null.
                  Read the rates as a description of the list, and open a gene to judge it.
                </p>
              </div>
            ) : null}

            {result.found.length ? (
              <div className="panel">
                <div className="panel-head">
                  <h2 className="panel-title">Genes, surest direction first</h2>
                  <button type="button" className="button-quiet" onClick={download}>
                    Download CSV
                  </button>
                </div>
                <div className="table-scroll">
                  <table className="list-table">
                    <thead>
                      <tr>
                        <th>Gene</th>
                        <th>Verdict</th>
                        <th className="num">g</th>
                        <th className="num">95% CI</th>
                        <th className="num" title="95% prediction interval: where one more study would land">
                          next study
                        </th>
                        <th className="num" title="independent GEO Series pooled; contrasts from one Series share subjects">series</th>
                        <th className="num" title="local false sign rate: the probability the direction is wrong">
                          P(wrong sign)
                        </th>
                        <th title="GenAge evidence code, where the gene is a curated human aging gene">
                          GenAge
                        </th>
                      </tr>
                    </thead>
                    <tbody>
                      {result.found.map((row) => (
                        <tr key={row.symbol}>
                          <td>
                            <button
                              type="button"
                              className="link-button mono"
                              onClick={() => onPick(row.symbol)}
                            >
                              {row.symbol}
                            </button>
                          </td>
                          <td className={row.verdict === "no_evidence" ? "null" : ""}>
                            <i className="dot" data-verdict={row.verdict} />
                            {VERDICT_LABEL[row.verdict]}
                          </td>
                          <td className="num">{fmtEffect(row.g)}</td>
                          <td className="num">{fmtInterval(row.ci_low, row.ci_high)}</td>
                          <td className="num">{fmtInterval(row.pi_low, row.pi_high)}</td>
                          <td className="num">{row.k}</td>
                          <td className="num">{fmtLfsr(row.lfsr)}</td>
                          <td className="mono">{evidence.get(row.symbol) ?? ""}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            ) : null}

            {result.missing.length ? (
              <div className="panel">
                <h2 className="panel-title">Not measured here ({fmtInt(result.missing.length)})</h2>
                <p className="mono missing-list">{result.missing.join(" ")}</p>
                <p className="note" style={{ maxWidth: "74ch" }}>
                  Not measured is not the same as no effect: these match no {species} gene with
                  three or more contrasts. Some may be retired names —{" "}
                  <a
                    className="structure-link"
                    href="https://www.genenames.org/tools/multi-symbol-checker/"
                    target="_blank"
                    rel="noreferrer"
                  >
                    check human symbols at HGNC ↗
                  </a>
                  .
                  {result.dates.length
                    ? ` ${fmtInt(result.dates.length)} ${result.dates.length === 1 ? "looks like a date" : "look like dates"}: spreadsheet software turns symbols such as SEPT1 and MARCH1 into "1-Sep" and "1-Mar". Export the list as plain text and paste it again.`
                    : ""}
                </p>
              </div>
            ) : null}
          </>
        ) : null}
      </div>
    </div>
  );
}
