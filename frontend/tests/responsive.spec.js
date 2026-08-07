import { expect, test } from "@playwright/test";

/**
 * Runs only under the `mobile` project (Pixel 5, 393px). "No mobile or
 * narrow-viewport work" was the first gap listed in the handoff, and the fix
 * was verified by looking at it once. A page that scrolls sideways is the one
 * responsive failure a reader notices immediately, so it is the one asserted.
 */
const VIEWS = ["gene", "landscape", "certainty", "panel", "mortality"];

for (const view of VIEWS) {
  test(`the ${view} view does not scroll the page sideways at 393px`, async ({ page }) => {
    await page.goto(`/?view=${view}`);
    await expect(
      page.locator(".chart svg, table tbody tr, .state-title").first(),
    ).toBeVisible({ timeout: 30_000 });

    const overflow = await page.evaluate(() => {
      const root = document.documentElement;
      return root.scrollWidth - root.clientWidth;
    });
    // Wide content is allowed to scroll inside .table-scroll / .chart; the
    // document body is not.
    expect(overflow).toBeLessThanOrEqual(1);
  });
}

test("wide tables scroll inside their own container", async ({ page }) => {
  await page.goto("/?view=panel");
  const scroller = page.locator(".table-scroll").first();
  await expect(scroller).toBeVisible({ timeout: 30_000 });

  const scrolls = await scroller.evaluate((el) => el.scrollWidth > el.clientWidth);
  const overflowX = await scroller.evaluate((el) => getComputedStyle(el).overflowX);
  expect(scrolls ? overflowX : "auto").toMatch(/auto|scroll/);
});
