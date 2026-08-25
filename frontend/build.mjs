// Bundles src/app.js (real ES module, imports Mol* directly) into app.js —
// the plain script the static frontend serves. Also copies Mol*'s
// precompiled stylesheet next to it. Run via `npm run build` (or
// `npm run watch` while developing).
import * as esbuild from "esbuild";
import { copyFileSync, mkdirSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const here = dirname(fileURLToPath(import.meta.url));
const watch = process.argv.includes("--watch");

mkdirSync(join(here, "vendor"), { recursive: true });
copyFileSync(
  join(here, "node_modules/molstar/build/viewer/molstar.css"),
  join(here, "vendor/molstar.css"),
);

const opts = {
  entryPoints: [join(here, "src/app.js")],
  outfile: join(here, "app.js"),
  bundle: true,
  format: "iife",
  platform: "browser",
  target: "es2020",
  sourcemap: true,
  logLevel: "info",
  define: { "process.env.NODE_ENV": '"production"' },
};

if (watch) {
  const ctx = await esbuild.context(opts);
  await ctx.watch();
  console.log("watching src/app.js for changes...");
} else {
  await esbuild.build(opts);
}
