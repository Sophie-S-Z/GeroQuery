"""Offline, batch, versioned ETL that regenerates the harmonized data layer.

Never on the request path. `build_fixtures.build_all()` writes the bundled,
redistributable demonstration slice; the production ETL swaps the synthetic
generators for live GEO/GTEx/recount3 harvests behind the same output schema.
"""
