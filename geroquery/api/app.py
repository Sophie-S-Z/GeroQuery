"""M7 api — FastAPI service layer.

Thin HTTP translation over :class:`GeroService`. Owns: routing, request
validation (Pydantic), the consistent error envelope, pagination, the
``format=json|csv|parquet`` switch, and OpenAPI docs. All domain logic lives
below this layer.
"""

from __future__ import annotations

import io

import pandas as pd
from fastapi import FastAPI, Query, Request, Response
from fastapi.responses import JSONResponse, StreamingResponse

from .. import __version__
from ..exceptions import GeroQueryError
from .schemas import (
    ClockApplyRequest,
    ClockCompareRequest,
    GeneSetRequest,
    ResilienceCSDRequest,
    ResilienceRecoveryRequest,
)
from .service import GeroService, get_service

API_PREFIX = "/v1"


def _matrix_from_request(service: GeroService, dataset_id, matrix, feature_cols=None):
    """Resolve a clock/resilience input to a DataFrame from either a dataset ID
    or an inline matrix payload."""
    if matrix is not None:
        return matrix.to_frame()
    if dataset_id is not None:
        df = service.store.get_dataset(dataset_id)
        return df[feature_cols] if feature_cols else df
    from ..exceptions import ClockInputError

    raise ClockInputError("Provide either 'dataset_id' or 'matrix'.")


def _formatted(rows: list[dict], fmt: str, filename: str) -> Response:
    """Return a list of dict rows as JSON, CSV, or Parquet."""
    fmt = (fmt or "json").lower()
    if fmt == "json":
        return JSONResponse(rows)
    df = pd.DataFrame(rows)
    if fmt == "csv":
        buf = io.StringIO()
        df.to_csv(buf, index=False)
        return StreamingResponse(
            iter([buf.getvalue()]),
            media_type="text/csv",
            headers={"Content-Disposition": f"attachment; filename={filename}.csv"},
        )
    if fmt == "parquet":
        bbuf = io.BytesIO()
        df.to_parquet(bbuf, index=False)
        bbuf.seek(0)
        return StreamingResponse(
            bbuf,
            media_type="application/octet-stream",
            headers={"Content-Disposition": f"attachment; filename={filename}.parquet"},
        )
    from ..exceptions import GeroQueryError as _E

    raise _E(f"Unsupported format {fmt!r}. Use json, csv, or parquet.")


def _paginate(items: list, limit: int | None, offset: int) -> list:
    if offset:
        items = items[offset:]
    if limit is not None:
        items = items[:limit]
    return items


def create_app(service: GeroService | None = None) -> FastAPI:
    app = FastAPI(
        title="GeroQuery API",
        version=__version__,
        description=(
            "Harmonized, multi-omic, cross-species, clock- and resilience-aware aging data."
        ),
    )

    # Lazily bind a service so importing the app never triggers a store build
    # until the first request (keeps import cheap and test setup explicit).
    _state: dict[str, GeroService | None] = {"service": service}

    def svc() -> GeroService:
        service_ = _state["service"]
        if service_ is None:
            service_ = get_service()
            _state["service"] = service_
        return service_

    @app.exception_handler(GeroQueryError)
    async def _domain_error(_: Request, exc: GeroQueryError):
        return JSONResponse(
            status_code=exc.http_status,
            content={"error": {"code": exc.code, "message": exc.message, "detail": exc.detail}},
        )

    # ---- health / version ------------------------------------------------

    @app.get("/healthz")
    def healthz():
        return {"status": "ok"}

    @app.get(f"{API_PREFIX}/version")
    def version():
        return svc().version()

    # ---- gene queries ----------------------------------------------------

    @app.get(f"{API_PREFIX}/gene/{{gene_id}}/signature")
    def gene_signature(
        gene_id: str,
        species: str | None = None,
        tissue: str | None = None,
        omic_layer: str | None = None,
        sex: str | None = None,
        limit: int | None = Query(None, ge=1),
        offset: int = Query(0, ge=0),
        format: str = Query("json", pattern="^(json|csv|parquet)$"),
    ):
        result = svc().gene_signature(gene_id, species, tissue, omic_layer, sex)
        result["signatures"] = _paginate(result["signatures"], limit, offset)
        if format != "json":
            return _formatted(result["signatures"], format, f"{gene_id}_signatures")
        return result

    @app.get(f"{API_PREFIX}/gene/{{gene_id}}/card")
    def gene_card(gene_id: str, species: str | None = None):
        return svc().gene_card(gene_id, species)

    @app.post(f"{API_PREFIX}/geneset/signature")
    def geneset_signature(req: GeneSetRequest):
        return svc().geneset_signature(req.genes, req.species, req.omic_layer)

    # ---- clocks ----------------------------------------------------------

    @app.get(f"{API_PREFIX}/clocks")
    def clocks():
        return {"clocks": svc().list_clocks()}

    @app.post(f"{API_PREFIX}/clock/apply")
    def clock_apply(req: ClockApplyRequest):
        df = _matrix_from_request(svc(), req.dataset_id, req.matrix)
        return svc().apply_clock(req.clock_id, df, req.chronological_age)

    @app.post(f"{API_PREFIX}/clock/compare")
    def clock_compare(req: ClockCompareRequest):
        df = _matrix_from_request(svc(), req.dataset_id, req.matrix)
        return {"results": svc().compare_clocks(req.clock_ids, df, req.chronological_age)}

    # ---- interventions ---------------------------------------------------

    @app.get(f"{API_PREFIX}/intervention/{{name}}")
    def intervention(name: str):
        return svc().intervention(name)

    # ---- resilience ------------------------------------------------------

    @app.post(f"{API_PREFIX}/resilience/csd")
    def resilience_csd(req: ResilienceCSDRequest):
        data = req.data.to_frame() if req.data is not None else None
        return svc().resilience_csd(
            dataset_id=req.dataset_id,
            data=data,
            biomarker_cols=req.biomarker_cols,
            age_col=req.age_col,
            n_strata=req.n_strata,
            longitudinal=req.longitudinal,
        )

    @app.post(f"{API_PREFIX}/resilience/recovery")
    def resilience_recovery(req: ResilienceRecoveryRequest):
        return svc().resilience_recovery(req.series)

    # ---- provenance ------------------------------------------------------

    @app.get(f"{API_PREFIX}/studies")
    def studies(
        limit: int | None = Query(None, ge=1),
        offset: int = Query(0, ge=0),
        format: str = Query("json", pattern="^(json|csv|parquet)$"),
    ):
        rows = _paginate(svc().studies(), limit, offset)
        if format != "json":
            return _formatted(rows, format, "studies")
        return {"studies": rows, "n": len(rows)}

    @app.get(f"{API_PREFIX}/sources")
    def sources():
        return {"sources": svc().sources()}

    @app.get(f"{API_PREFIX}/datasets")
    def datasets():
        return {"datasets": svc().datasets()}

    return app


# Module-level ASGI app for `uvicorn geroquery.api.app:app`.
app = create_app()
