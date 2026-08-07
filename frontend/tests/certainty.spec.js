import { expect, test } from "@playwright/test";
import { pooledRows, watchForErrors } from "./helpers.js";

/**
 * The certainty view, and the data contract behind it.
 *
 * The chart tests assert the view renders and reads correctly. The data-contract
 * test is the one that earns its place: two of this repo's silent UI failures
 * were defects in the exported Parquet that rendered without complaint, and bug
 * #20 was a *third* — a verdict computed at a precision the reader never saw.
 * Its first fix covered four of six interval columns and missed the two added in
 * the same commit, so the assertions below are written over the column list
 * rather than over columns named by hand.
 */

const VIEW = "/?view=certainty";

/** Every bound the export ships. Enumerated once; see the note above. */
const BOUND_COLUMNS = ["ci_low", "ci_high", "ci_low_dl", "ci_high_dl", "pi_low", "pi_high"];

test.describe("certainty view", () => {
  test("draws every gene against its false sign rate", async ({ page }) => {
    const problems = watchForErrors(page);
    await page.goto(VIEW);

    const chart = page.locator(".chart").first();
    await expect(chart.locator("svg")).toBeVisible({ timeout: 30_000 });

    // One mark per pooled gene, plus the emphasised and labelled ones. If this
    // collapses to a handful, the join or the filter is wrong and the chart
    // still looks like a chart.
    const dots = chart.locator("svg circle");
    expect(await dots.count()).toBeGreaterThan(10_000);
    expect(problems).toEqual([]);
  });

  test("states the gap between intervals excluding zero and confident directions", async ({
    page,
  }) => {
    await page.goto(VIEW);
    // `.certainty-summary`, not `.note` — the summary renders only once the
    // data lands, so `.note` first-match resolves to the table caption below
    // and the assertions run against the wrong element.
    const note = page.locator(".certainty-summary");
    await expect(note).toBeVisible({ timeout: 30_000 });

    // The whole argument of the view: strictly more genes have an interval
    // excluding zero than have a small lfsr. Read the numbers off the page and
    // assert the inequality, rather than asserting the sentence exists.
    // Thousands-separated integers only, so "0.05" cannot be mistaken for one.
    const [excludes, total, confident] = [...(await note.textContent()).matchAll(/\b(\d{1,3}(?:,\d{3})+|\d{3,})\b/g)]
      .slice(0, 3)
      .map((m) => Number(m[1].replace(/,/g, "")));

    expect(confident).toBeGreaterThan(0);
    expect(total).toBeGreaterThan(excludes);
    expect(excludes).toBeGreaterThan(confident);
  });

  test("a labelled gene navigates to its own readout", async ({ page }) => {
    await page.goto(VIEW);
    const first = page.locator("table tbody tr td button").first();
    await expect(first).toBeVisible({ timeout: 30_000 });
    const symbol = (await first.textContent()).trim();

    await first.click();
    await expect(page.locator("h1.gene-symbol")).toHaveText(symbol, { timeout: 30_000 });
  });

  test("the landscape table is ranked by the shrunken effect, not the raw one", async ({
    page,
  }) => {
    // The winner's curse fix, asserted where it matters. Mouse is the species
    // where the two orderings diverge hardest: its former top row by raw g had
    // a 46% posterior chance of the wrong sign.
    await page.goto("/?view=landscape&species=mouse");
    const rows = page.locator("table tbody tr");
    await expect(rows.first()).toBeVisible({ timeout: 30_000 });

    const shrunk = await page.locator("table tbody tr td:nth-child(3)").allTextContents();
    const values = shrunk.map((cell) => Math.abs(Number(cell)));
    expect(values.length).toBeGreaterThan(5);

    // Descending in the shrunken column is the claim. Ranked by raw g this
    // column would be visibly out of order — that is exactly what the bug was.
    const sorted = [...values].sort((a, b) => b - a);
    expect(values).toEqual(sorted);
  });
});

test.describe("pooled data contract", () => {
  test("no interval bound ships as negative zero", async () => {
    // Bug #20's residual. `round(-1e-9, 4)` is `-0.0`, which prints as "-0.000"
    // and compares as *not less than zero* in JavaScript, so a row can disagree
    // with itself across the language boundary. Two genes shipped `pi_high`
    // that way after the first fix covered only the four confidence bounds.
    const rows = await pooledRows();
    expect(rows.length).toBeGreaterThan(0);

    for (const column of BOUND_COLUMNS) {
      const negativeZero = rows.filter(
        (row) => row[column] === 0 && Object.is(row[column], -0),
      );
      expect(negativeZero, `${column} carries negative zero`).toEqual([]);
    }
  });

  test("every gene carries a posterior, and it is a probability", async () => {
    const rows = await pooledRows();

    const missing = rows.filter((row) => row.lfsr == null || row.g_shrunk == null);
    expect(missing).toEqual([]);

    const outOfRange = rows.filter((row) => row.lfsr < 0 || row.lfsr > 1);
    expect(outOfRange).toEqual([]);
  });

  test("shrinkage pulls toward zero and never past it", async () => {
    // The defining property of a posterior mean under a prior that is unimodal
    // at zero. A sign flip or an overshoot means the conjugate update is wrong,
    // and neither would look wrong plotted.
    const rows = await pooledRows();

    const overshot = rows.filter((row) => Math.abs(row.g_shrunk) > Math.abs(row.g) + 1e-6);
    expect(overshot).toEqual([]);

    const flipped = rows.filter(
      (row) => row.g_shrunk !== 0 && Math.sign(row.g_shrunk) !== Math.sign(row.g),
    );
    expect(flipped).toEqual([]);
  });
});
