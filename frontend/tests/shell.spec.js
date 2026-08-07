import { expect, test } from "@playwright/test";
import { ACCESSION, watchForErrors } from "./helpers.js";

const VIEWS = ["gene", "landscape", "panel", "mortality"];

/**
 * Bugs #9 and #10 were whole tabs that rendered but could never load anything —
 * a Streamlit branch on a key that was always null, and a radio whose options no
 * longer matched its conditions. Both left a page that returned HTTP 200 and
 * showed nothing. So the shell test is: every view must put real content on the
 * page, and must do it without the console complaining.
 */
test.describe("shell", () => {
  for (const view of VIEWS) {
    test(`the ${view} view renders content and logs nothing`, async ({ page }) => {
      const problems = watchForErrors(page);
      await page.goto(`/?view=${view}`);

      await expect(page.locator(`.tab[aria-selected="true"]`)).toHaveText(new RegExp(view, "i"));
      // Something that took data to draw, not just a heading the bundle carries.
      await expect(
        page.locator(".chart svg, table tbody tr, .state-title").first(),
      ).toBeVisible({ timeout: 30_000 });
      await expect(page.locator(".skeleton")).toHaveCount(0);
      expect(problems).toEqual([]);
    });
  }

  test("the masthead counts come from the data, not from the source", async ({ page }) => {
    await page.goto("/");
    const stats = page.locator(".masthead-stats");
    await expect(stats).toBeVisible();
    await expect(stats).toContainText(/\d{1,3}(,\d{3})+ effect sizes/);
    await expect(stats).toContainText(/\d+ contrasts \/ \d+ series/);

    const meta = await page.evaluate(() => fetch("/data/meta.json").then((r) => r.json()));
    await expect(stats).toContainText(meta.n_genes.toLocaleString("en-US"));
    await expect(page.locator(".mode-flag")).toHaveAttribute("data-mode", meta.signature_mode);
  });

  test("switching theme redraws the charts rather than leaving stale colours", async ({ page }) => {
    await page.goto("/?view=landscape");
    const chart = page.locator(".chart svg").first();
    await expect(chart).toBeVisible({ timeout: 30_000 });

    const before = await chart.getAttribute("style");
    await page.locator(".theme-toggle").click();
    await expect(page.locator(".chart svg").first()).toBeVisible();
    const after = await page.locator(".chart svg").first().getAttribute("style");
    expect(after).not.toBe(before);
  });

  test("the landscape table lists symbols, never Ensembl accessions", async ({ page }) => {
    await page.goto("/?view=landscape&species=human&dir=up");
    const genes = page.locator("table tbody tr td:first-child button");
    await expect(genes.first()).toBeVisible({ timeout: 30_000 });

    const symbols = await genes.allTextContents();
    expect(symbols.length).toBeGreaterThan(5);
    // Bug #17 rendered 94.4% of genes as their own accession. One is too many.
    expect(symbols.filter((s) => ACCESSION.test(s))).toEqual([]);
  });

  test("the view tabs are reachable and operable from the keyboard", async ({ page }) => {
    await page.goto("/");
    const landscape = page.locator(".tab", { hasText: "Landscape" });
    await landscape.focus();
    await expect(landscape).toBeFocused();
    await page.keyboard.press("Enter");
    await expect(page).toHaveURL(/view=landscape/);
  });
});
