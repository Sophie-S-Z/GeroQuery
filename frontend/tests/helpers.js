import { readFile } from "node:fs/promises";
import { parquetReadObjects } from "hyparquet";
import { expect } from "@playwright/test";

/** An Ensembl accession rendered where a gene symbol belongs. Bug #17. */
export const ACCESSION = /^ENS[A-Z]*G?\d{6,}/;

/**
 * Open a view and wait for its data, collecting anything the console complained
 * about on the way.
 *
 * Returns the collected problems rather than asserting on them, so a test can
 * decide whether a given failure is the thing it is about. `shell.spec.js`
 * asserts the list is empty for every view; the others use it as context when
 * they fail.
 */
export function watchForErrors(page) {
  const problems = [];
  page.on("console", (message) => {
    if (message.type() === "error") problems.push(`console: ${message.text()}`);
  });
  page.on("pageerror", (error) => problems.push(`pageerror: ${error.message}`));
  page.on("requestfailed", (request) =>
    problems.push(`requestfailed: ${request.url()} ${request.failure()?.errorText}`),
  );
  return problems;
}

/** Wait until the gene readout has resolved to either a result or a stated absence. */
export async function waitForGene(page) {
  await expect(page.locator("h1.gene-symbol, .state-title").first()).toBeVisible({
    timeout: 30_000,
  });
}

/** Every text label on a Plot y-axis, in draw order. */
export async function axisLabels(chart, axis = "y") {
  return chart.locator(`svg [aria-label="${axis}-axis tick label"] text`).allTextContents();
}

/**
 * The deployed pooled corpus, read in Node straight off `dist`.
 *
 * Not through the page: two of the silent failures this suite exists for were
 * defects in the exported data that rendered without complaint, so the file is
 * worth asserting on as a file.
 */
export async function pooledRows() {
  const file = await readFile(new URL("../dist/data/pooled.parquet", import.meta.url));
  const rows = await parquetReadObjects({
    file: file.buffer.slice(file.byteOffset, file.byteOffset + file.byteLength),
    utf8: true,
  });
  return rows;
}
