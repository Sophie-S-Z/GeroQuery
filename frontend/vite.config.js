import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// No dev proxy and no API base: this app has no backend. It reads two Parquet
// files from /data with hyparquet, which is why it can be a static bundle on a
// CDN and cannot go down or exceed a free tier.
export default defineConfig({
  plugins: [react()],
  // GitHub Pages serves project sites from /<repo>/, so the workflow sets
  // BASE_PATH=/GeroQuery/. Cloudflare Pages serves from the root and needs
  // nothing. `DATA_BASE` in src/lib/db.js reads import.meta.env.BASE_URL, so
  // the Parquet files resolve under either.
  base: process.env.BASE_PATH || "/",
  build: {
    outDir: "dist",
    target: "es2022",
  },
  server: { port: 5173 },
});
