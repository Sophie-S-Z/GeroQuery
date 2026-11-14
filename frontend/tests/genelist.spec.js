import { readFile } from "node:fs/promises";
import { expect, test } from "@playwright/test";
import { watchForErrors } from "./helpers.js";

async function curatedJson() {
  return JSON.parse(await readFile(new URL("../dist/data/curated.json", import.meta.url), "utf8"));
}

const pct = (x) => `${Math.round(x * 100)}%`;

test.describe("gene list", () => {
  /**
   * The browser counts a pasted list; Python counted GenAge when it built the
   * data. Loading GenAge into the list must reproduce the Python numbers
   * exactly, or the two paths have drifted.
   */
  test("loading GenAge reproduces the counts the pipeline shipped", async ({ page }) => {
    const problems = watchForErrors(page);
    const curated = await curatedJson();
    expect(curated.available).toBe(true);

    await page.goto("/?view=list");
    await page.getByRole("button", { name: /^Try GenAge/ }).click();

    await expect(page.locator("h1.section-title")).toHaveText(
      `${curated.n_measured} of ${curated.n_listed} genes measured`,
      { timeout: 30_000 },
    );
    const resolved = curated.verdicts.increases + curated.verdicts.decreases;
    const lede = page.locator(".list-lede");
    await expect(lede).toContainText(`${resolved} have a pooled direction with age (${pct(curated.rate_resolved)})`);
    await expect(lede).toContainText(`against ${pct(curated.background.rate_resolved)}`);
    await expect(lede).toContainText(`${curated.n_pi_excludes_zero} clear the prediction interval`);
    await expect(page.locator(".list-table tbody tr")).toHaveCount(curated.n_measured);
    if (curated.not_measured.length) {
      await expect(page.locator(".missing-list")).toHaveText(curated.not_measured.join(" "));
    } else {
      await expect(page.locator(".missing-list")).toHaveCount(0);
    }
    expect(page.url()).toContain("genes=");
    expect(problems).toEqual([]);
  });

  test("a pasted list names what it could not match, including spreadsheet dates", async ({
    page,
  }) => {
    await page.goto("/?view=list");
    await page.locator("#gene-list").fill("FOXO3, cdkn2a\nNOTAGENE  1-Sep\nFOXO3");
    await page.getByRole("button", { name: "Check list" }).click();

    await expect(page.locator("h1.section-title")).toHaveText("2 of 4 genes measured", { timeout: 30_000 });
    await expect(page.locator(".missing-list")).toHaveText("NOTAGENE 1-Sep");
    await expect(page.locator(".panel", { hasText: "Not measured" })).toContainText("1 looks like a date");
    await expect(page.locator(".list-table tbody tr")).toHaveCount(2);

    // The link reproduces the check, and a gene opens its own page.
    await page.reload();
    await expect(page.locator("h1.section-title")).toHaveText("2 of 4 genes measured", { timeout: 30_000 });
    await page.locator(".list-table tbody button", { hasText: "CDKN2A" }).click();
    await expect(page.locator("h1.gene-symbol")).toHaveText("CDKN2A", { timeout: 30_000 });
  });

  test("the CSV holds exactly the rows on screen", async ({ page }) => {
    await page.goto("/?view=list&genes=FOXO3,CDKN2A,IGF1");
    await expect(page.locator(".list-table tbody tr")).toHaveCount(3, { timeout: 30_000 });

    const [download] = await Promise.all([
      page.waitForEvent("download"),
      page.getByRole("button", { name: "Download CSV" }).click(),
    ]);
    expect(download.suggestedFilename()).toMatch(/^geroquery-gene-list-human-.+\.csv$/);
    const lines = (await readFile(await download.path(), "utf8")).trim().split("\n");
    expect(lines[0]).toBe("symbol,species,verdict,g,ci_low,ci_high,pi_low,pi_high,k,n_contrasts,lfsr,genage_evidence");
    expect(lines).toHaveLength(4);
  });

  test("the overview's curated-gene finding opens the GenAge list", async ({ page }) => {
    const curated = await curatedJson();
    await page.goto("/");
    await page.locator(".findings button", { hasText: "curated aging genes" }).click();
    await expect(page).toHaveURL(/view=list/);
    await expect(page.locator("h1.section-title")).toHaveText(
      `${curated.n_measured} of ${curated.n_listed} genes measured`,
      { timeout: 30_000 },
    );
  });
});

test.describe("gene page actions", () => {
  test("the citation carries the estimate as shown and the data version", async ({
    page,
    context,
  }) => {
    await context.grantPermissions(["clipboard-read", "clipboard-write"]);
    await page.goto("/?view=gene&gene=CDKN2A&species=human");
    const pooled = page.locator(".pooled-line");
    await expect(pooled).toBeVisible({ timeout: 30_000 });

    await page.getByRole("button", { name: "Copy citation" }).click();
    await expect(page.getByRole("button", { name: "Citation copied" })).toBeVisible();
    const text = await page.evaluate(() => navigator.clipboard.readText());
    const meta = await page.evaluate(() => fetch("/data/meta.json").then((r) => r.json()));
    const g = (await pooled.textContent()).match(/g = ([+-]\d\.\d{3})/)[1];

    expect(text).toContain(`CDKN2A (human): pooled g = ${g}`);
    expect(text).toContain(`data version ${meta.data_version}`);
    expect(text).toContain("gene=CDKN2A");
  });

  test("GenAge is linked under the name it lists a renamed gene by", async ({ page }) => {
    const curated = await curatedJson();
    const [old, current] = Object.entries(curated.renamed)[0];
    await page.goto(`/?view=gene&gene=${current}&species=human`);
    const link = page.getByRole("link", { name: /GenAge/ });
    await expect(link).toHaveAttribute("href", `https://genomics.senescence.info/genes/entry.php?hgnc=${old}`, {
      timeout: 30_000,
    });
    await expect(page.getByRole("link", { name: /GTEx/ })).toBeVisible();
  });

  test("a mouse gene gets no human-only links", async ({ page }) => {
    await page.goto("/?view=gene&gene=CDKN2A&species=mouse");
    await expect(page.getByRole("link", { name: /NCBI Gene/ })).toHaveAttribute("href", /Mus%20musculus/, {
      timeout: 30_000,
    });
    await expect(page.getByRole("link", { name: /GTEx|GenAge/ })).toHaveCount(0);
  });

  test("Ctrl+K jumps to gene search from another view", async ({ page }) => {
    await page.goto("/?view=panel");
    await expect(page.locator(".panel").first()).toBeVisible({ timeout: 30_000 });
    await page.keyboard.press("Control+k");
    await expect(page.locator("#gene-search")).toBeFocused();
  });
});
