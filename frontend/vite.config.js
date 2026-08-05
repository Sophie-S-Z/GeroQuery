import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// No dev proxy and no API base: this app has no backend. It reads two Parquet
// files from /data with hyparquet, which is why it can be a static bundle on a
// CDN and cannot go down or exceed a free tier.
export default defineConfig({
  plugins: [react()],
  build: {
    outDir: "dist",
    target: "es2022",
  },
  server: { port: 5173 },
});
