import { readFile } from "node:fs/promises";
import { expect, test } from "@playwright/test";
import { axisLabels, watchForErrors } from "./helpers.js";

async function influenceJson() {
  const raw = await readFile(new URL("../dist/data/influence_corpus.json", import.meta.url), "utf8");
  return JSON.parse(raw);
}

test.describe("corpus influence (§4.5)", () => {
  test("ranks the contrasts by how hard they pull, and accounts for every one", async ({
    page,
  }) => {
    const problems = watchForErrors(page);
    await page.goto("/?view=panel");
    await expect(page.locator(".panel-title", { hasText: "bends the corpus" })).toBeVisible({
      timeout: 30_000,
    });

    const ranking = await influenceJson();
    expect(ranking.length).toBeGreaterThan(1);

    // The head of each species' ranking is drawn (at most 20 bars), and each
    // caption states the full count so the cut is visible, never silent.
    const studies = JSON.parse(
      await readFile(new URL("../dist/data/studies.json", import.meta.url), "utf8"),
    );
    const species = Object.fromEntries(studies.map((s) => [s.study_id, s.species]));
    const perSpecies = {};
    for (const r of ranking) perSpecies[species[r.study_id]] = (perSpecies[species[r.study_id]] ?? 0) + 1;
    const panel = page.locator(".panel:has(h2:text-is('Which contrast bends the corpus'))");
    const labels = await axisLabels(panel.locator(".chart"));
    expect(labels.length).toBe(
      Object.values(perSpecies).reduce((sum, n) => sum + Math.min(20, n), 0),
    );
    for (const [name, n] of Object.entries(perSpecies)) {
      await expect(panel.locator(".facet-title", { hasText: name })).toContainText(String(n));
    }

    // Ranked by the mean shift per gene measured, not the total, which mostly
    // ranks array size (coverage runs 1,007 to 20,897 genes).
    const shifts = ranking.map((r) => r.mean_abs_shift);
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
    await expect(page.locator("th", { hasText: "mean shift" })).toBeVisible();
    // The table prints the top contrast's mean shift to three decimals.
    await expect(page.locator("tbody")).toContainText(top.mean_abs_shift.toFixed(3));
  });
});
