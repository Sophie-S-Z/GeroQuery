# Deploying GeroQuery

The front end is a static bundle. There is no server, no database, and no
environment variable to set — which is why it can be hosted free, permanently,
without a credit card, and cannot go down or exceed a quota.

`frontend/dist/` is 8.2 MB total. Largest single file is 5.1 MB, well under the
25 MiB per-file limit of every host below.

---

## Fastest path — Cloudflare Pages via Wrangler (~2 minutes)

Recommended. Unlimited bandwidth, never sleeps, free, no card.

Run these yourself — `wrangler login` opens a browser window for the OAuth
consent, which is why it cannot be done for you.

```bash
cd frontend
npm install
npm run build

npx wrangler login                       # opens your browser, approve once
npx wrangler pages project create geroquery --production-branch main
npx wrangler pages deploy dist --project-name geroquery
```

The last command prints the live URL, something like
`https://geroquery.pages.dev`. That is it — no further configuration.

To redeploy after changing anything:

```bash
make frontend                            # re-export data + rebuild
cd frontend && npx wrangler pages deploy dist --project-name geroquery
```

> In Claude Code you can run any of these inline by prefixing with `!`, e.g.
> `! npx wrangler login` — the output lands in the conversation.

---

## Alternative — Cloudflare Pages via the dashboard (git-connected)

Better if you want every push to `main` to redeploy automatically.

1. Push the repo to GitHub (it already points at
   `github.com/Sophie-S-Z/GeroQuery`).
2. Go to **Cloudflare dashboard → Workers & Pages → Create → Pages → Connect to Git**.
3. Pick the `GeroQuery` repository.
4. Set:

   | Setting | Value |
   |---|---|
   | Production branch | `main` |
   | Framework preset | None |
   | Build command | `npm run build` |
   | Build output directory | `dist` |
   | Root directory | `frontend` |

5. **Settings → Environment variables → Production**, add `NODE_VERSION` = `20`.
6. Save and Deploy.

`frontend/public/_headers` is picked up automatically and sets the caching and
`Accept-Ranges` headers. Nothing else needs configuring.

**Why the data has to be committed.** Cloudflare's build container runs
`npm run build` only — it has no Python, no store, and no id maps, so it cannot
run `python -m geroquery.etl.build_frontend_data`. `frontend/public/data/` is
therefore committed (6.8 MB) as the deployable artifact. `meta.json` carries the
data version that identifies exactly which corpus a deployed site is showing.

---

## Alternative — GitHub Pages (no new account)

Free, and the repo is already on GitHub. Slightly worse than Cloudflare
(bandwidth is soft-capped at 100 GB/month and it is a little slower), but it
needs no third-party signup.

A workflow is committed at `.github/workflows/deploy-pages.yml`. To turn it on:

1. **Repo → Settings → Pages → Build and deployment → Source: GitHub Actions.**
2. Push to `main`, or run the workflow manually from the Actions tab.

The site lands at `https://sophie-s-z.github.io/GeroQuery/`.

The workflow sets `BASE_PATH=/GeroQuery/` because GitHub Pages serves project
sites from a subpath. `vite.config.js` reads it, and `DATA_BASE` in
`src/lib/db.js` already honours `import.meta.env.BASE_URL`, so the Parquet files
resolve correctly under the subpath. **Cloudflare needs no such setting** — it
serves from the root.

---

## Local preview

```bash
make frontend           # export data + build
cd frontend
npm run preview         # http://localhost:4173
```

Or serve the built output with anything static:

```bash
npx serve -s frontend/dist -l 4173
```

Use the built output rather than `npm run dev` when checking anything that
matters. The dev server transforms modules on the fly and does not exercise the
asset hashing, the `_headers` file, or the production bundle — and this repo has
a documented history of UI bugs that only appeared in the real artifact.

---

## Rebuilding the data

The site reads five files from `frontend/public/data/`. Regenerate them with:

```bash
make frontend-data
```

That requires a built store. Full sequence from a clean checkout:

```bash
pip install -e ".[dev,ui]"
make data              # ~679 MB, ~30 min, needs network
make crosslayer        # the mortality cohort
make frontend-data
```

**`make idmap` must have run at least once** or 94% of genes render as raw
Ensembl accessions instead of symbols — the exporter inverts the bulk id maps
that step writes. It raises if no map exists at all, but it cannot tell a stale
map from a complete one.

Without `make data`, the store builds from the committed curated slice and the
site shows 40,585 rows instead of 485,905. That is not silent: `meta.json`
records `signature_mode`, and the masthead displays a **CURATED SLICE** flag in
place of **FULL CORPUS**.

---

## What is deployed where

| Component | Where it runs | Cost |
|---|---|---|
| Front end | Cloudflare Pages / GitHub Pages, static | free |
| Signature + mortality data | Same CDN, as Parquet/JSON | free |
| Living-evidence loop | GitHub Actions, monthly cron | free on public repos |
| MCP server | Locally, via `geroquery-mcp` | n/a |
| FastAPI service | Not deployed; optional | — |

The API is deliberately not part of the deployment. Everything the site needs is
precomputed, so hosting it would add a running service, a free-tier ceiling, and
a cold start for no gain. If you later want it live — for the clock library on
uploaded cohorts, which genuinely needs Python — Google Cloud Run's always-free
tier (2M requests/month, scales to zero) is the option that survives; see
`docs/ROADMAP.md` §3.2 for the comparison and why Hugging Face Spaces no longer
qualifies.

---

## Connecting the MCP server

Not a deployment — it runs locally and is registered with your client.

```bash
pip install -e ".[mcp]"
```

Then add to your MCP client config (for Claude Code, `~/.claude.json`):

```json
{
  "mcpServers": {
    "geroquery": {
      "command": "geroquery-mcp"
    }
  }
}
```

Five tools become available: `gene_aging_signature`, `geneset_aging_signature`,
`intervention_effect`, `list_studies`, `data_provenance`. Every response carries
the confidence interval and contrast count in the primary payload, and returns
`no_evidence` rather than a point estimate when the interval crosses zero.
