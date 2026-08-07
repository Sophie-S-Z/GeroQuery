import { defineConfig, devices } from "@playwright/test";

/**
 * Playwright runs against the **built** `dist`, never the dev server.
 *
 * Every silent UI failure this repo has shipped (#9, #10, #17, #18, #19) was
 * invisible to a green Python suite because nothing rendered the page. Two of
 * them — the raw-accession labels and the merged forest rows — depended on the
 * exported Parquet, which the dev server and the build read identically but
 * which no unit test touches at all. So the harness has to be the real artifact:
 * `vite preview` over `dist`, the same bytes Cloudflare serves.
 *
 * Chromium only. A second engine would double CI time to re-check rendering
 * behaviour that is not where the bugs have been; the bugs have been in what
 * the data says, and that is engine-independent.
 */
export default defineConfig({
  testDir: "./tests",
  fullyParallel: true,
  forbidOnly: !!process.env.CI,
  retries: process.env.CI ? 1 : 0,
  reporter: process.env.CI ? [["github"], ["list"]] : [["list"]],
  use: {
    baseURL: "http://127.0.0.1:4173",
    trace: "retain-on-failure",
  },
  projects: [
    {
      name: "desktop",
      use: { ...devices["Desktop Chrome"], viewport: { width: 1440, height: 900 } },
      // The responsive spec asserts a 393px-specific property; running it at
      // 1440 would report a pass that means nothing.
      testIgnore: /responsive\.spec\.js/,
    },
    { name: "mobile", use: { ...devices["Pixel 5"] }, testMatch: /responsive\.spec\.js/ },
  ],
  webServer: {
    command: "npm run preview",
    url: "http://127.0.0.1:4173",
    reuseExistingServer: !process.env.CI,
    // The corpus is 7 MB of Parquet; the first build after a data rebuild is
    // the slow part, not the server.
    timeout: 120_000,
  },
});
