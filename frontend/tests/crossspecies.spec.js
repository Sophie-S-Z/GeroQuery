import { readFile } from "node:fs/promises";
import { expect, test } from "@playwright/test";
import { watchForErrors } from "./helpers.js";

async function crossSpeciesJson() {
  const raw = await readFile(new URL("../dist/data/crossspecies.json", import.meta.url), "utf8");
  return JSON.parse(raw);
}

test.describe("cross-species (§4.2)", () => {
  test("plots every ortholog pair and shows the Deming fit, not an OLS one", async ({ page }) => {
    const problems = watchForErrors(page);
    await page.goto("/?view=crossspecies");
    await expect(page.locator("h1.section-title")).toHaveText(/Human vs mouse/, { timeout: 30_000 });

    const payload = await crossSpeciesJson();
    expect(payload.available).toBe(true);

    // One dot per pair (Plot renders dots as <circle>), plus a visible chart.
    const chart = page.locator(".chart svg").first();
    await expect(chart).toBeVisible();
    const circles = await chart.locator("circle").count();
    expect(circles).toBe(payload.n_pairs);

    // The honest headline is on the page: the resolved count and zero discordant.
    const resolved = payload.quadrants.both_up + payload.quadrants.both_down + payload.quadrants.discordant;
    await expect(page.locator(".panel", { hasText: "What the picture says" })).toContainText(
      String(resolved),
    );
    expect(problems).toEqual([]);
  });

  test("the resolved-genes table matches the shipped verdicts", async ({ page }) => {
    await page.goto("/?view=crossspecies");
    await expect(page.locator("h1.section-title")).toBeVisible({ timeout: 30_000 });

    const payload = await crossSpeciesJson();
    const resolved = payload.pairs.filter(
      (p) => p.verdict_human !== "no_evidence" && p.verdict_mouse !== "no_evidence",
    );
    const table = page.locator(".panel:has(h2:text-is('Genes resolved in both species')) tbody tr");
    await expect(table).toHaveCount(resolved.length);
  });
});
