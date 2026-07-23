import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// Dev server proxies /v1 to the FastAPI backend so the showcase and API share an
// origin during development. In production set VITE_API_BASE to the API URL.
export default defineConfig({
  plugins: [react()],
  server: {
    proxy: {
      "/v1": "http://localhost:8000",
    },
  },
});
