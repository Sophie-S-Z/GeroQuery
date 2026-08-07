/**
 * Data layer — two Parquet files, read in the browser, no backend.
 *
 * The first cut of this used DuckDB-WASM so that queries could range-request
 * individual row groups. That was the wrong trade and Cloudflare Pages' 25 MiB
 * per-file limit is what made it obvious: the DuckDB wasm binaries are 34 MB
 * and 39 MB, about seven times the size of the 6.6 MB corpus they existed to
 * query. Range requests only pay when the data dwarfs the engine.
 *
 * So: hyparquet, a ~15 kB pure-JS Parquet reader. Both files are fetched whole
 * exactly once, lazily, and every query after that is a local array operation —
 * no network, no worker, no wasm. Total transfer is 6.8 MB, under a fifth of
 * the old engine alone, and the site still has no server.
 *
 * The hard rule this module holds: **no statistics are computed here.** Pooled
 * effects, intervals, verdicts and hazard ratios are all computed in Python by
 * the same estimators the API serves; this file filters, sorts and bins rows
 * that already carry their answers. A DerSimonian-Laird re-implementation in
 * JavaScript would be a second estimator for one quantity, and two estimators
 * drift. Filtering is rendering; pooling is deriving.
 */

import { parquetReadObjects } from "hyparquet";

/** Base URL for the data directory, honouring a non-root deploy. */
const DATA_BASE = `${import.meta.env.BASE_URL ?? "/"}data`.replace(/\/{2,}/g, "/");

const cache = new Map();

/**
 * Parquet int64 columns (here: `k`, the contrast count) decode to BigInt.
 *
 * BigInt looks harmless until it reaches something that is not expecting it:
 * `JSON.stringify` throws outright, mixing it with a Number in arithmetic
 * throws a TypeError, and Plot's scales silently produce NaN. Normalizing once
 * at the boundary is the only place this is cheap — every call site downstream
 * would otherwise have to remember, and one that forgot would fail at render
 * time rather than at load time.
 *
 * These counts are all far below Number.MAX_SAFE_INTEGER, so the conversion is
 * lossless.
 */
function normalizeBigInts(row) {
  let converted = null;
  for (const key in row) {
    if (typeof row[key] === "bigint") {
      converted ??= { ...row };
      converted[key] = Number(row[key]);
    }
  }
  return converted ?? row;
}

/**
 * Load and memoize one Parquet file as plain row objects.
 *
 * The promise is cached, not the result, so concurrent callers during the
 * initial load share one fetch instead of racing three of them.
 */
function table(name) {
  if (!cache.has(name)) {
    cache.set(
      name,
      (async () => {
        const url = new URL(`${DATA_BASE}/${name}`, window.location.href).href;
        // One plain fetch, not hyparquet's asyncBufferFromUrl.
        //
        // The ranged reader issues a HEAD plus a widening series of Range GETs
        // per file, and because every row group is read anyway those ranges
        // overlap: measured at 8 requests and 7.58 MB transferred for 6.8 MB of
        // files. Fetching the whole buffer once is fewer requests, less bytes,
        // and drops the dependency on the host honouring Range at all.
        const response = await fetch(url);
        if (!response.ok) throw new Error(`HTTP ${response.status}`);
        const file = await response.arrayBuffer();
        const rows = await parquetReadObjects({ file, utf8: true });
        return rows.map(normalizeBigInts);
      })().catch((error) => {
        // Do not memoize a failure: a transient network error should not
        // poison every later query for the life of the page.
        cache.delete(name);
        throw new Error(`Could not load ${name}: ${error.message}`);
      }),
    );
  }
  return cache.get(name);
}

const pooledTable = () => table("pooled.parquet");
const contrastsTable = () => table("contrasts.parquet");

/** Static JSON fetched whole — all of these are a few KB. */
export async function json(name) {
  const response = await fetch(`${DATA_BASE}/${name}`);
  if (!response.ok) throw new Error(`Could not load ${name} (${response.status})`);
  return response.json();
}

// --- queries -----------------------------------------------------------------

/**
 * Symbol prefix search, for the query rail.
 *
 * Shorter symbols first: typing "TP5" should surface TP53 before TP53BP2,
 * because the shorter match is nearly always the one being looked for.
 */
export async function searchGenes(term, limit = 40) {
  const needle = term.toUpperCase();
  const rows = await pooledTable();
  const hits = rows.filter((row) => row.symbol?.toUpperCase().startsWith(needle));
  hits.sort(
    (a, b) =>
      a.symbol.length - b.symbol.length ||
      a.symbol.localeCompare(b.symbol) ||
      a.species.localeCompare(b.species),
  );
  return hits.slice(0, limit);
}

/** Every pooled row for one symbol — one per species. */
export async function pooledForSymbol(symbol) {
  const needle = symbol.toUpperCase();
  const rows = await pooledTable();
  return rows
    .filter((row) => row.symbol?.toUpperCase() === needle)
    .sort((a, b) => a.species.localeCompare(b.species));
}

/** Per-study contrasts behind a pooled estimate. This is the forest plot. */
export async function contrastsForSymbol(symbol, species) {
  const needle = symbol.toUpperCase();
  const rows = await contrastsTable();
  return rows
    .filter((row) => row.symbol?.toUpperCase() === needle && row.species === species)
    .sort((a, b) => a.effect_size - b.effect_size);
}

/**
 * Histogram of pooled effects across the whole panel, split by verdict.
 *
 * Binned here rather than in the chart so the mark gets a few hundred rows
 * instead of 41,983 — this is the one view whose job is to show the shape of
 * the entire corpus at once, and Plot should not be asked to aggregate it.
 */
export async function effectDistribution(species, binWidth = 0.1) {
  const rows = await pooledTable();
  const counts = new Map();
  for (const row of rows) {
    if (row.species !== species) continue;
    if (row.g < -3 || row.g > 3) continue;
    const bin = Math.round(Math.floor(row.g / binWidth) * binWidth * 1e4) / 1e4;
    const key = `${row.verdict}|${bin}`;
    const entry = counts.get(key);
    if (entry) entry.n += 1;
    else counts.set(key, { verdict: row.verdict, bin, n: 1 });
  }
  return [...counts.values()].sort((a, b) => a.bin - b.bin);
}

/**
 * The strongest replicated effects, for the landscape table.
 *
 * `rank` picks which estimate orders the table, and the choice is not cosmetic:
 *
 * `"shrunk"` (default)
 *     The corpus posterior mean. Ranking by the raw pooled effect selects on
 *     significance and then sorts on an extreme of a noisy statistic, which is
 *     the winner's curse — the top of that table is enriched for exactly the
 *     genes that will shrink on replication. It is not hypothetical here: the
 *     mouse table's former top row, `Kif21a`, has an interval excluding zero
 *     and a **46% posterior probability that its direction is wrong**.
 * `"raw"`
 *     The old ordering, kept so the two can be compared rather than asserted.
 */
export async function topGenes({ species, direction, minK = 3, limit = 60, rank = "shrunk" }) {
  const rows = await pooledTable();
  const excludesZero =
    direction === "down" ? (row) => row.ci_high < 0 : (row) => row.ci_low > 0;
  // Fall back to the raw effect for any row with no posterior — a species with
  // too few genes to fit a corpus prior still has to sort.
  const key = (row) => (rank === "raw" || row.g_shrunk == null ? row.g : row.g_shrunk);
  return rows
    .filter((row) => row.species === species && row.k >= minK && excludesZero(row))
    .sort((a, b) => (direction === "down" ? key(a) - key(b) : key(b) - key(a)))
    .slice(0, limit);
}

/**
 * Every pooled estimate that carries a posterior, for the volcano.
 *
 * Returns raw rows: the chart bins and draws them, and every value it reads —
 * `g`, `lfsr`, `verdict` — was computed in Python. No thresholding happens
 * here, deliberately. Selecting features by effect size *and* significance is
 * the double filter that breaks the false-discovery-rate guarantee, so this
 * view offers no "large and significant" selector to break it with.
 */
export async function volcanoPoints(species) {
  const rows = await pooledTable();
  return rows.filter((row) => row.species === species && row.lfsr != null);
}

/** Corpus-level counts, computed rather than typed into the masthead. */
export async function panelSummary() {
  const rows = await pooledTable();
  const bySpecies = new Map();
  for (const row of rows) {
    let entry = bySpecies.get(row.species);
    if (!entry) {
      entry = { species: row.species, n_pooled: 0, n_up: 0, n_down: 0, n_null: 0, widths: [] };
      bySpecies.set(row.species, entry);
    }
    entry.n_pooled += 1;
    if (row.verdict === "increases") entry.n_up += 1;
    else if (row.verdict === "decreases") entry.n_down += 1;
    else entry.n_null += 1;
    entry.widths.push(row.ci_high - row.ci_low);
  }
  return [...bySpecies.values()]
    .map(({ widths, ...entry }) => {
      widths.sort((a, b) => a - b);
      const mid = Math.floor(widths.length / 2);
      const median =
        widths.length % 2 ? widths[mid] : (widths[mid - 1] + widths[mid]) / 2;
      return { ...entry, median_ci_width: Math.round(median * 1000) / 1000 };
    })
    .sort((a, b) => a.species.localeCompare(b.species));
}
