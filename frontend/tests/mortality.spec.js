import { expect, test } from "@playwright/test";
import { axisLabels, watchForErrors } from "./helpers.js";

test.describe("mortality result", () => {
  test("every clock gets its own row, in the chart and in the table", async ({ page }) => {
    const problems = watchForErrors(page);
    await page.goto("/?view=mortality");

    const chart = page.locator(".chart").first();
    await expect(chart.locator("svg")).toBeVisible({ timeout: 30_000 });

    const labels = await axisLabels(chart);
    expect(new Set(labels).size).toBe(labels.length);

    const rows = page.locator("table tbody tr");
    const clocks = await rows.locator("td:first-child").allTextContents();
    expect(clocks.length).toBe(labels.length);

    /**
     * Bug #19: `/Age$|Mort$|AgeMort$/` matches leftmost-first, so `GrimAgeMort`
     * became "Grim" while `GrimAge2Mort` became "GrimAge2" — the same clock
     * family labelled at two different levels of abbreviation.
     *
     * Uniqueness does not catch that; both labels are distinct. What catches it
     * is the contract: a chart label is its clock name with at most the literal
     * suffix "Mort" removed, and nothing else. Under the bug the leftover for
     * GrimAgeMort is "AgeMort" and for PhenoAge it is "Age".
     */
    labels.forEach((label, index) => {
      const clock = clocks[index];
      expect(clock.startsWith(label)).toBe(true);
      expect(clock.slice(label.length)).toMatch(/^(|Mort)$/);
    });
    expect(problems).toEqual([]);
  });

  test("the cohort description matches the exported payload", async ({ page }) => {
    await page.goto("/?view=mortality");
    await expect(page.locator(".panel-title").first()).toBeVisible({ timeout: 30_000 });

    const data = await page.evaluate(() =>
      fetch("/data/crosslayer.json").then((r) => r.json()),
    );
    test.skip(!data.available, "cross-layer cohort not built in this checkout");

    const note = page.locator(".panel .note").first();
    await expect(note).toContainText(data.n_subjects.toLocaleString("en-US"));
    await expect(note).toContainText(data.n_deaths.toLocaleString("en-US"));
    await expect(note).toContainText(String(data.followup_years_median));
  });

  test("hazard ratios are shown with their intervals, never bare", async ({ page }) => {
    await page.goto("/?view=mortality");
    const rows = page.locator("table tbody tr");
    await expect(rows.first()).toBeVisible({ timeout: 30_000 });

    for (const row of await rows.all()) {
      const cells = await row.locator("td").allTextContents();
      // HR, then a bracketed interval beside it. The product's whole claim is
      // that a point estimate without an interval is not an answer.
      expect(cells[2]).toMatch(/^\d+\.\d+$/);
      expect(cells[3]).toMatch(/\[.*,.*\]/);
    }
  });

  test("the reported hazard ratios are the weighted ones, and say so", async ({ page }) => {
    await page.goto("/?view=mortality");
    await expect(page.locator("table tbody tr").first()).toBeVisible({ timeout: 30_000 });

    const data = await page.evaluate(() =>
      fetch("/data/crosslayer.json").then((r) => r.json()),
    );
    test.skip(!data.available, "cross-layer cohort not built in this checkout");
    expect(data.survey?.weighted).toBe(true);

    // The population claim has to be on the page, not only in the payload —
    // an unlabelled survey-weighted number reads as a sample number.
    await expect(page.locator(".panel .note").first()).toContainText(
      data.survey.population_size.toLocaleString("en-US"),
    );

    // Column 3 is the weighted HR and column 5 the unweighted one; a build that
    // wired them the other way round would still render plausibly.
    const first = page.locator("table tbody tr").first();
    const cells = await first.locator("td").allTextContents();
    expect(Number(cells[2])).toBeCloseTo(data.clocks[0].hr, 3);
    expect(Number(cells[4])).toBeCloseTo(data.clocks[0].hr_sample, 3);
    expect(data.clocks[0].hr).toBeGreaterThan(data.clocks[0].hr_sample);
  });

  test("the limitations are rendered, not filed away", async ({ page }) => {
    await page.goto("/?view=mortality");
    const caveats = page.locator("ul.caveats li");
    await expect(caveats.first()).toBeVisible({ timeout: 30_000 });
    expect(await caveats.count()).toBeGreaterThan(3);
  });
});
