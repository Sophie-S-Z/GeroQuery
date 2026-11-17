import { expect, test } from "@playwright/test";

/**
 * The perturbation layer must read as a separate, caveated surface. These tests
 * pin the parts that stop it being mistaken for an aging verdict.
 */

test("the perturbation view states that neither cell line is aged tissue", async ({ page }) => {
  await page.goto("/?view=perturbation");
  await expect(page.getByText(/Neither cell line is aged tissue/i)).toBeVisible();
  await expect(page.getByRole("heading", { name: /two cell lines mostly disagree/i })).toBeVisible();
});

test("coverage names the genes no screen has ever switched off", async ({ page }) => {
  await page.goto("/?view=perturbation");
  const row = page.getByRole("row", { name: /Absent from both screens/i });
  await expect(row).toBeVisible();
  // A real count, not a placeholder.
  await expect(row.locator("td.num")).toHaveText(/^[\d,]+$/);
});

test("replicating genes are reachable and carry a direction", async ({ page }) => {
  await page.goto("/?view=perturbation");
  const table = page.locator(".panel", { hasText: "replicates across both cell lines" });
  await expect(table.locator("tbody tr").first()).toBeVisible();
  await expect(table.getByText(/switching it off looks (more|less) like age/).first()).toBeVisible();
});

test("a gene never switched off says so rather than showing a score", async ({ page }) => {
  await page.goto("/?view=gene&gene=CDKN1A&species=human");
  const panel = page.locator(".panel", { hasText: "Switching this gene off" });
  await expect(panel).toBeVisible();
  await expect(panel.getByText(/has not been switched off in either screen/i)).toBeVisible();
});

test("a tested gene shows both cell lines with the caveat attached", async ({ page }) => {
  await page.goto("/?view=gene&gene=DICER1&species=human");
  const panel = page.locator(".panel", { hasText: "Switching this gene off" });
  await expect(panel.getByRole("row", { name: /RPE1/ })).toBeVisible();
  await expect(panel.getByRole("row", { name: /K562/ })).toBeVisible();
  await expect(panel.getByText(/lead for an experiment|means little|disagree/i)).toBeVisible();
});

test("the mouse gene view offers no knockdown panel, since the screens are human", async ({
  page,
}) => {
  await page.goto("/?view=gene&gene=Dicer1&species=mouse");
  await expect(page.locator(".panel", { hasText: "Switching this gene off" })).toHaveCount(0);
});
