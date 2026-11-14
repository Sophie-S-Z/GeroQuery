import { readFile } from "node:fs/promises";
import { expect, test } from "@playwright/test";
import { pooledRows, waitForGene } from "./helpers.js";

/**
 * Microarray and RNA-seq contrasts are pooled apart in Python; the page shows
 * those two estimates as shipped and counts nothing itself.
 */
test.describe("cross-platform replication", () => {
  test("a gene pooled on both platforms shows each platform's estimate", async ({ page }) => {
    const rows = await pooledRows();
    const gene = rows
      .filter((r) => r.species === "human" && r.k_array != null && r.k_rnaseq != null)
      .sort((a, b) => Number(b.k) - Number(a.k))[0];
    test.skip(!gene, "no gene is pooled on both platforms in this build");

    await page.goto(`/?view=gene&gene=${gene.symbol}&species=human`);
    await waitForGene(page);
    const table = page.locator(".panel", { hasText: "pooled apart" }).locator("tbody tr");
    await expect(table).toHaveCount(2);
    const fmt = (x) => `${x >= 0 ? "+" : "−"}${Math.abs(x).toFixed(2)}`;
    // Parquet int64 arrives in Node as BigInt; the page normalizes it.
    await expect(table.nth(0)).toContainText(String(gene.k_array));
    await expect(table.nth(1)).toContainText(String(gene.k_rnaseq));
    // d3-format writes the minus sign as U+2212.
    await expect(table.nth(0).locator("td").nth(1)).toHaveText(fmt(gene.g_array));
    await expect(table.nth(1).locator("td").nth(1)).toHaveText(fmt(gene.g_rnaseq));
  });

  test("the replication finding opens the replicated genes", async ({ page }) => {
    const meta = JSON.parse(await readFile(new URL("../dist/data/meta.json", import.meta.url), "utf8"));
    test.skip(!meta.platforms?.n_replicated, "nothing replicates in this build");
    const [species, genes] = Object.entries(meta.platforms.replicated).sort(
      (a, b) => b[1].length - a[1].length,
    )[0];

    await page.goto("/");
    const finding = page.locator(".findings button", { hasText: "replicate across arrays" });
    await expect(finding).toContainText(meta.platforms.n_replicated.toLocaleString("en-US"));
    await finding.click();
    await expect(page.locator("h1.section-title")).toHaveText(
      `${genes.length.toLocaleString("en-US")} of ${genes.length.toLocaleString("en-US")} genes measured`,
      { timeout: 30_000 },
    );
    await expect(page.locator('.species-switch button[aria-pressed="true"]')).toHaveText(species);
  });
});

test("the readout counts independent Series and contrasts, as shipped", async ({ page }) => {
  const rows = await pooledRows();
  const gene = rows.find(
    (r) => r.symbol === "CDKN2A" && r.species === "human" && r.n_contrasts != null,
  );
  test.skip(!gene, "export predates series-level pooling");
  await page.goto("/?view=gene&gene=CDKN2A&species=human");
  await waitForGene(page);
  await expect(page.locator(".pooled-line")).toContainText(
    `k = ${Number(gene.k)} series (${Number(gene.n_contrasts)} contrasts)`,
  );
});
