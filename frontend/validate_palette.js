/**
 * Validate GeroQuery's diverging polarity palette — in both themes, every
 * theme block, against the actual shipped `src/styles.css`.
 *
 * `VISUAL_AND_STATISTICAL_PLAN.md` §1 rule 4 says colour is *validated, not
 * chosen*, and §6 is the audit this script performs. This is NOT the generic
 * categorical validator: GeroQuery's palette is a three-way **diverging** scale
 * — `--up` (warm), `--down` (cool), `--null` (a deliberately warm neutral grey)
 * — and the categorical lightness-band / chroma-floor checks reject a correct
 * diverging midpoint *by design* (§6 records this so nobody "fixes" the neutral
 * by saturating it). So the checks here are the ones that matter for polarity:
 *
 *   1. Pairwise separation of the three verdict fills, in OKLab ΔE×100, under
 *      normal vision and simulated protan/deutan/tritan (Machado 2009, sev 1.0).
 *   2. Contrast of each fill against its own surface (WCAG ≥ 3:1).
 *
 * The gate reflects what §6 actually adjudicated. The confusable pair — the one
 * whose collapse "destroys the distinction the whole tool exists to draw" — is
 * **null ↔ down** ("interval crosses zero" vs "interval below zero"); the two
 * poles **up ↔ down** are the other hard pair. Both are hard-gated at the ΔE 15
 * normal-vision floor and a CVD floor of 8. The third pair, **up ↔ null**, is a
 * warm pole against a warm neutral and is geometrically the closest of the
 * three — a diverging midpoint cannot sit far from *both* poles at once. It is
 * reported and WARN-gated: it must clear the CVD floor, and below the normal
 * floor it is legal only because §5 makes secondary encoding (direct labels,
 * legend, position on the zero-aligned forest plot) mandatory for every verdict,
 * so identity is never colour-alone.
 *
 * Reads every theme block in styles.css: the light `:root`, the explicit
 * `:root[data-theme="dark"]`, and the `@media (prefers-color-scheme: dark)`
 * block — because a half-applied fix to one dark path and not the other is
 * exactly the §6 defect, and only checking every block catches it.
 *
 * Usage:  node frontend/validate_palette.js [path/to/styles.css]
 * Exit 0 if every theme passes; 1 on any hard FAIL; 2 on a usage/parse error.
 */

import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";

// -- thresholds (OKLab ΔE ×100; WCAG for contrast) -----------------------------
const NORMAL_FLOOR = 15.0; // hard floor for the audited pairs, unsimulated vision
const CVD_FLOOR = 8.0; // hard floor for the audited pairs, min(protan, deutan)
const CONTRAST_MIN = 3.0; // each fill vs its surface

// Machado, Oliveira & Fernandes (2009) CVD transforms at severity 1.0 (linear RGB).
const MACHADO = {
  protan: [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
  deutan: [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
  tritan: [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]],
};

// -- colour conversions (sRGB → linear → OKLab; WCAG relative luminance) --------
const hex2srgb = (h) => { h = h.trim().replace(/^#/, ""); return [0, 2, 4].map((i) => parseInt(h.slice(i, i + 2), 16) / 255); };
const s2lin = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
const lin = (h) => hex2srgb(h).map(s2lin);
const relLum = (h) => { const [r, g, b] = lin(h); return 0.2126 * r + 0.7152 * g + 0.0722 * b; };
const contrast = (a, b) => { const [hi, lo] = [relLum(a), relLum(b)].sort((x, y) => y - x); return (hi + 0.05) / (lo + 0.05); };

function oklabFromLin([r, g, b]) {
  const l = Math.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b);
  const m = Math.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b);
  const s = Math.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b);
  return [
    0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
    1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
    0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s,
  ];
}
function simulate(h, kind) {
  const [r, g, b] = lin(h), M = MACHADO[kind], c = (x) => Math.max(0, Math.min(1, x));
  return [c(M[0][0] * r + M[0][1] * g + M[0][2] * b), c(M[1][0] * r + M[1][1] * g + M[1][2] * b), c(M[2][0] * r + M[2][1] * g + M[2][2] * b)];
}
function deltaE(h1, h2, kind) {
  const a = oklabFromLin(kind ? simulate(h1, kind) : lin(h1));
  const b = oklabFromLin(kind ? simulate(h2, kind) : lin(h2));
  return 100 * Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

// -- extract the theme token blocks from styles.css ----------------------------
// Each block is a `{ ... }` body; we read --up/--down/--null and the surface
// (--paper) out of whichever blocks define the polarity tokens.
function themesFromCss(css) {
  const token = (body, name) => {
    const m = body.match(new RegExp(`--${name}\\s*:\\s*(#[0-9a-fA-F]{6})`));
    return m ? m[1] : null;
  };
  // Pull out brace-balanced bodies for each selector we care about, in order.
  const blocks = [
    { label: "light (:root)", selector: /:root\s*\{/g, mode: "light" },
    { label: 'dark (:root[data-theme="dark"])', selector: /:root\[data-theme="dark"\]\s*\{/g, mode: "dark" },
    { label: "dark (@media prefers-color-scheme)", selector: /:root:not\(\[data-theme="light"\]\)\s*\{/g, mode: "dark" },
  ];
  const themes = [];
  for (const { label, selector, mode } of blocks) {
    const m = selector.exec(css);
    if (!m) continue;
    // Walk braces from the opening `{` to its match.
    let depth = 0, start = m.index + m[0].length - 1, i = start;
    for (; i < css.length; i++) {
      if (css[i] === "{") depth++;
      else if (css[i] === "}" && --depth === 0) break;
    }
    const body = css.slice(start + 1, i);
    const up = token(body, "up"), down = token(body, "down"), nul = token(body, "null"), surface = token(body, "paper");
    if (up && down && nul && surface) themes.push({ label, mode, up, down, null: nul, surface });
  }
  return themes;
}

// -- the audit -----------------------------------------------------------------
const PAIRS = [
  ["down", "null", "hard", "confusable verdicts — the pair §6 exists to protect"],
  ["up", "down", "hard", "the two poles"],
  ["up", "null", "warn", "warm pole vs warm neutral — relies on mandatory secondary encoding (§5)"],
];

function auditTheme(t) {
  const rows = [];
  let hardFail = false;
  for (const [a, b, gate, why] of PAIRS) {
    const normal = deltaE(t[a], t[b]);
    const cvd = Math.min(deltaE(t[a], t[b], "protan"), deltaE(t[a], t[b], "deutan"));
    const tritan = deltaE(t[a], t[b], "tritan");
    let state;
    if (gate === "hard") {
      state = normal >= NORMAL_FLOOR && cvd >= CVD_FLOOR ? "PASS" : "FAIL";
      if (state === "FAIL") hardFail = true;
    } else {
      state = normal >= NORMAL_FLOOR && cvd >= CVD_FLOOR ? "PASS" : cvd >= CVD_FLOOR ? "WARN" : "FAIL";
      if (state === "FAIL") hardFail = true;
    }
    rows.push({ pair: `${a}↔${b}`, state, normal: normal.toFixed(1), cvdWorst: cvd.toFixed(1), tritan: tritan.toFixed(1), why });
  }
  const contrasts = [];
  for (const k of ["up", "down", "null"]) {
    const cr = contrast(t[k], t.surface);
    const state = cr >= CONTRAST_MIN ? "PASS" : "FAIL";
    if (state === "FAIL") hardFail = true;
    contrasts.push({ fill: k, ratio: `${cr.toFixed(2)}:1`, state });
  }
  return { rows, contrasts, ok: !hardFail };
}

export { deltaE, contrast, themesFromCss, auditTheme };

// -- CLI -----------------------------------------------------------------------
function isMain() {
  return process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url);
}

if (isMain()) {
  const here = dirname(fileURLToPath(import.meta.url));
  const path = process.argv[2] ? resolve(process.argv[2]) : resolve(here, "src", "styles.css");
  let css;
  try {
    css = readFileSync(path, "utf8");
  } catch (e) {
    console.error(`cannot read ${path}: ${e.message}`);
    process.exit(2);
  }
  const themes = themesFromCss(css);
  if (!themes.length) {
    console.error(`no theme token blocks with --up/--down/--null/--paper found in ${path}`);
    process.exit(2);
  }
  console.log(`\nGeroQuery palette audit — ${path}`);
  console.log(`floors: normal ΔE ≥ ${NORMAL_FLOOR}, CVD ΔE ≥ ${CVD_FLOOR}, contrast ≥ ${CONTRAST_MIN}:1  (OKLab ΔE×100, Machado 2009)\n`);
  let allOk = true;
  for (const t of themes) {
    const { rows, contrasts, ok } = auditTheme(t);
    allOk = allOk && ok;
    console.log(`  ${t.label}  up ${t.up} · down ${t.down} · null ${t.null} · surface ${t.surface}`);
    for (const r of rows) {
      console.log(`    [${r.state.padEnd(4)}] ${r.pair.padEnd(10)} normal ${r.normal.padStart(5)}  cvd ${r.cvdWorst.padStart(5)}  tritan ${r.tritan.padStart(5)}  — ${r.why}`);
    }
    for (const c of contrasts) {
      console.log(`    [${c.state.padEnd(4)}] contrast ${c.fill.padEnd(5)} ${c.ratio} vs surface`);
    }
    console.log("");
  }
  console.log(allOk ? "→ ALL THEMES PASS\n" : "→ FAILED — fix the marked checks\n");
  process.exit(allOk ? 0 : 1);
}
