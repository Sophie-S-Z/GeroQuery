import { readFile } from "node:fs/promises";
import { expect, test } from "@playwright/test";
import { watchForErrors } from "./helpers.js";

async function metaJson() {
  const raw = await readFile(new URL("../dist/data/meta.json", import.meta.url), "utf8");
  return JSON.parse(raw);
}

test.describe("overview (landing)", () => {
  test("is the default view and leads with the prediction-interval headline", async ({ page }) => {
    const problems = watchForErrors(page);
    await page.goto("/");

    // Default view, no query string needed.
    await expect(page.locator('.tab[aria-selected="true"]')).toHaveText(/overview/i);

    const meta = await metaJson();
    const pi = meta.prediction_interval;
    // The honest headline number is the PI count, and it is the one in the hero.
    await expect(page.locator(".overview-head")).toContainText(
      pi.n_pi_excludes_zero.toLocaleString("en-US"),
    );
    // Both counts are stated, PI smaller than CI, so the gap is legible.
    await expect(page.locator(".overview-lede")).toContainText(
      pi.n_ci_excludes_zero.toLocaleString("en-US"),
    );
    expect(pi.n_pi_excludes_zero).toBeLessThan(pi.n_ci_excludes_zero);

    // The funnel chart draws (its own data-driven svg).
    await expect(page.locator(".overview-hero .chart svg")).toBeVisible();
    expect(problems).toEqual([]);
  });

  test("a finding row navigates into its view", async ({ page }) => {
    await page.goto("/");
    const row = page.locator(".findings button", { hasText: "pulls hardest on the evidence" });
    await expect(row).toBeVisible({ timeout: 30_000 });
    await row.click();
    await expect(page).toHaveURL(/view=panel/);
    await expect(page.locator(".panel-title", { hasText: "bends the corpus" })).toBeVisible();
  });

  test("the wordmark returns home from any view", async ({ page }) => {
    await page.goto("/?view=gene");
    await page.locator(".wordmark").click();
    await expect(page).toHaveURL(/^[^?]*\/?$|view=overview/);
    await expect(page.locator(".overview-head")).toBeVisible({ timeout: 30_000 });
  });

  test("the search box is the primary action and opens the gene", async ({ page }) => {
    await page.goto("/");
    const field = page.getByLabel("Look up a gene");
    await field.fill("FOXO3");
    await field.press("Enter");
    await expect(page).toHaveURL(/view=gene/);
    await expect(page).toHaveURL(/gene=FOXO3/);
    await expect(page.locator("h1.gene-symbol")).toHaveText("FOXO3", { timeout: 30_000 });
  });

  test("/ focuses gene search from any view", async ({ page }) => {
    await page.goto("/?view=panel");
    await expect(page.locator(".tab[aria-selected=\"true\"]")).toHaveText(/panel/i);
    await page.keyboard.press("/");
    await expect(page).toHaveURL(/view=gene/);
    await expect(page.locator("#gene-search")).toBeFocused();
  });
});
