import { readFile } from "node:fs/promises";
import { expect, test } from "@playwright/test";
import { axisLabels, watchForErrors } from "./helpers.js";

async function influenceJson() {
  const raw = await readFile(new URL("../dist/data/influence_corpus.json", import.meta.url), "utf8");
  return JSON.parse(raw);
}

test.describe("corpus influence (§4.5)", () => {
  test("ranks every contrast by how much it moves the whole corpus", async ({ page }) => {
    const problems = watchForErrors(page);
    await page.goto("/?view=panel");
    await expect(page.locator(".panel-title", { hasText: "bends the corpus" })).toBeVisible({
      timeout: 30_000,
    });

    const ranking = await influenceJson();
    expect(ranking.length).toBeGreaterThan(1);

    // One bar per contrast, and the axis is ranked most-influential first.
    const chart = page.locator(".panel:has(h2:text-is('Which contrast bends the corpus')) .chart");
    const labels = await axisLabels(chart);
    expect(labels.length).toBe(ranking.length);

    // Ranked: the shipped file is sorted by total_abs_shift descending.
    const shifts = ranking.map((r) => r.total_abs_shift);
    for (let i = 1; i < shifts.length; i++) expect(shifts[i]).toBeLessThanOrEqual(shifts[i - 1]);

    expect(problems).toEqual([]);
  });

  test("the table carries the same corpus-shift numbers as the chart", async ({ page }) => {
    await page.goto("/?view=panel");
    await expect(page.locator(".panel-title", { hasText: "bends the corpus" })).toBeVisible({
      timeout: 30_000,
    });

    const ranking = await influenceJson();
    const top = ranking[0];
    // Every chart needs a table; the study table is it, carrying the shift column.
    await expect(page.locator("th", { hasText: "corpus shift" })).toBeVisible();
    // The table prints the top contrast's shift to one decimal.
    await expect(page.locator("tbody")).toContainText(top.total_abs_shift.toFixed(1));
  });
});
