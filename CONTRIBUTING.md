# Contributing

Thanks for your interest in GeroQuery.

## Setup

```bash
pip install -e ".[dev,ui]"
pre-commit install
make test
```

## Adding a data source

Write **one adapter** in `geroquery/sources/` implementing the `SourceAdapter`
interface (`capabilities()`, `license()`, and the relevant `fetch_*`). Nothing
else in the system should need to change. If the source is controlled-access,
set `cacheable=False` and `redistributable=False` on its capabilities/licence —
the store will refuse to cache it and the licence test will confirm.

## Conventions

- Deep modules: keep interfaces small and stable; hide complexity inside.
- Lower modules never import higher ones (`ui → api → services → sources/idmap`).
- Tests assert **external behavior against ground truth**, not internals.
- Scientific modules (M1, M3, M5, M6) get unit tests; plumbing (M2, M4, M7, M8)
  gets integration/contract tests.
- `ruff` + `black` + `mypy` must pass (`make ci`).

## Scientific correctness

If you touch M3 (harmonize), M5 (clocks), or M6 (resilience), add a test that
pins the new behavior to a known-truth case (a planted effect, a published age,
an engineered tipping point). Correctness here is the whole point.
