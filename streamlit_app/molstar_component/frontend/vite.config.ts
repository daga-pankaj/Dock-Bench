import { defineConfig } from "vite";

// Streamlit serves the built component from a sub-path, so assets must be
// referenced relatively.
export default defineConfig({
  base: "./",
  build: {
    outDir: "dist",
    chunkSizeWarningLimit: 12000,
  },
  define: { "process.env.NODE_ENV": JSON.stringify("production") },
  test: { environment: "node", include: ["src/**/*.test.ts"] },
} as any);
