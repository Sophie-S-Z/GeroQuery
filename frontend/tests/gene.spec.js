import { expect, test } from "@playwright/test";
import { ACCESSION, axisLabels, pooledRows, waitForGene, watchForErrors } from "./helpers.js";

test.describe("gene readout", () => {
  test("the default view resolves to a real gene with an interval", async ({ page }) => {
    const problems = watchForErrors(page);
    await page.goto("/");
    await waitForGene(page);

    await expect(page.locator("h1.gene-symbol")).toHaveText(/^[A-Z0-9-]+$/);
    // g = 0.072 [-0.203, +0.347] · k = 14 · p = ...
    await expect(page.locator(".pooled-line")).toContainText(/g = [-+0-9.]/);
    await expect(page.locator(".pooled-line")).toContainText(/k = \d+/);
    expect(problems).toEqual([]);
  });

  /**
   * Bug #18. One GDS can yield two contrasts, and the forest plot used to label
   * both with the bare accession — Plot's ordinal y-scale then merged them into
   * a single row, drawing two measurements as one.
   *
   * Mouse Cdkn2a is the fixture because GDS5226 is the only accession in the
   * panel that yields two contrasts (bone marrow and epididymal adipocyte) and
   * Cdkn2a is measured in both. Human CDKN2A, the obvious choice, touches no
   * duplicated accession at all — this test passed against the reintroduced bug
   * until the fixture was picked from the data rather than from memory.
   */
  test("the forest plot draws one row per contrast and never merges two", async ({ page }) => {
    await page.goto("/?gene=Cdkn2a&species=mouse&view=gene");
    await waitForGene(page);

    const chart = page.locator(".chart").first();
    await expect(chart.locator("svg")).toBeVisible();

    const labels = await axisLabels(chart);
    const contrastRows = await page
      .locator('.panel:has(h2:text-is("Contrast table")) tbody tr')
      .count();

    // Two rows from GDS5226, so the fixture can actually expose the merge.
    expect(contrastRows).toBeGreaterThan(1);
    expect(labels.filter((l) => l.startsWith("GDS5226")).length).toBe(2);
    expect(labels).toContain("POOLED");
    expect(new Set(labels).size).toBe(labels.length);
    // Exactly one row per contrast, plus the pooled summary. Not "at most".
    expect(labels.length).toBe(contrastRows + 1);
  });

  test("a gene that was never measured says so, rather than showing an empty chart", async ({
    page,
  }) => {
    await page.goto("/?gene=NOTAREALGENE&species=human&view=gene");
    await waitForGene(page);

    await expect(page.locator(".state-title")).toHaveText("Not measured");
    await expect(page.locator(".state")).toContainText("different from finding no effect");
    await expect(page.locator("h1.gene-symbol")).toHaveCount(0);
  });

  /** Bug #17: the search rail could not find 94% of the corpus by name. */
  test("search finds genes by symbol and never offers a raw accession", async ({ page }) => {
    await page.goto("/");
    await waitForGene(page);

    await page.locator("#gene-search").fill("CDKN");
    const results = page.locator(".result-symbol");
    await expect(results.first()).toBeVisible({ timeout: 15_000 });

    const symbols = await results.allTextContents();
    expect(symbols.length).toBeGreaterThan(1);
    for (const symbol of symbols) expect(symbol).not.toMatch(ACCESSION);
    expect(symbols.some((s) => s.startsWith("CDKN"))).toBe(true);
  });

  test("picking a result moves the URL, so a result is linkable", async ({ page }) => {
    await page.goto("/");
    await waitForGene(page);

    await page.locator("#gene-search").fill("TP53");
    const first = page.locator(".result").first();
    await expect(first).toBeVisible({ timeout: 15_000 });
    const picked = await first.locator(".result-symbol").textContent();
    await first.click();

    await expect(page).toHaveURL(new RegExp(`gene=${picked}`));
    await expect(page.locator("h1.gene-symbol, .state-title")).toContainText(
      new RegExp(`${picked}`),
    );
  });

  test("the species switch changes the pooled estimate, not just the label", async ({ page }) => {
    await page.goto("/?gene=CDKN1A&species=human&view=gene");
    await waitForGene(page);

    const switcher = page.locator(".species-switch button", { hasText: "mouse" });
    test.skip((await switcher.count()) === 0, "CDKN1A has no mouse pool in this build");

    const before = await page.locator(".pooled-line").textContent();
    await switcher.click();
    await expect(page).toHaveURL(/species=mouse/);
    await expect(page.locator(".gene-meta")).toContainText("mouse");
    await expect(page.locator(".pooled-line")).not.toHaveText(before ?? "");
  });
});

test.describe("intervals", () => {
  /**
   * The prediction interval is the one statistic on this page that answers a
   * different question from the one beside it, so a build that dropped it would
   * leave the page saying something narrower than it means.
   */
  test("a pooled gene shows both the confidence and the prediction interval", async ({ page }) => {
    await page.goto("/?gene=CDKN2A&species=human&view=gene");
    await waitForGene(page);

    await expect(page.locator(".pooled-line")).toContainText(/\[.*,.*\]/);
    const prediction = page.locator(".note", { hasText: "prediction interval" }).first();
    await expect(prediction).toBeVisible();
    await expect(prediction).toContainText(/Next study \[.*,.*\]/);
  });

  /**
   * A data-contract check on the deployed Parquet, run in Node rather than in
   * the page: two of this repo's silent UI failures were defects in the export
   * that rendered perfectly, so the file itself is worth asserting on directly.
   */
  test("the exported corpus carries honest intervals", async () => {
    const rows = await pooledRows();
    expect(rows.length).toBeGreaterThan(1000);

    // Never narrower than the interval every other tool reports.
    const narrower = rows.filter(
      (r) => r.ci_high - r.ci_low < r.ci_high_dl - r.ci_low_dl - 1e-9,
    );
    expect(narrower).toEqual([]);

    // A prediction interval contains the confidence interval it was built from.
    const withPi = rows.filter((r) => r.pi_low != null);
    expect(withPi.length).toBe(rows.length); // every pooled row has k >= 3
    const inverted = withPi.filter((r) => r.pi_high - r.pi_low < r.ci_high - r.ci_low);
    expect(inverted).toEqual([]);

    // The verdict is the reported interval's, not the DL one's.
    const mislabelled = rows.filter((r) => {
      const expected =
        r.ci_low > 0 ? "increases" : r.ci_high < 0 ? "decreases" : "no_evidence";
      return r.verdict !== expected;
    });
    expect(mislabelled).toEqual([]);
  });
});
