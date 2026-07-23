"""Pydantic request/response schemas — the HTTP contract for /v1.

Only request bodies and a couple of response wrappers need declaring; most
responses are assembled dicts from the service and documented via examples.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
from pydantic import BaseModel, Field, field_validator


class GeneSetRequest(BaseModel):
    genes: list[str] = Field(
        ..., min_length=1, description="Gene identifiers (symbol/Ensembl/Entrez)."
    )
    species: str | None = None
    omic_layer: str | None = None

    @field_validator("genes")
    @classmethod
    def _cap(cls, v):
        if len(v) > 500:
            raise ValueError("A gene set may contain at most 500 genes.")
        return v


class MatrixPayload(BaseModel):
    """A feature matrix as row records plus optional sample IDs.

    ``records`` is a list of {feature: value} dicts, one per sample — the natural
    JSON shape of a user-uploaded expression/methylation/clinical matrix.
    """

    records: list[dict[str, Any]] = Field(..., min_length=1)
    sample_ids: list[str] | None = None

    @field_validator("records")
    @classmethod
    def _cap_rows(cls, v):
        if len(v) > 20000:
            raise ValueError("Matrix exceeds 20000 rows; use the batch/ETL path instead.")
        return v

    def to_frame(self) -> pd.DataFrame:
        df = pd.DataFrame(self.records)
        if self.sample_ids is not None:
            if len(self.sample_ids) != len(df):
                raise ValueError("sample_ids length must match number of records.")
            df.index = self.sample_ids
        return df


class ClockApplyRequest(BaseModel):
    clock_id: str
    dataset_id: str | None = None
    matrix: MatrixPayload | None = None
    chronological_age: list[float] | None = None


class ClockCompareRequest(BaseModel):
    clock_ids: list[str] = Field(..., min_length=1)
    dataset_id: str | None = None
    matrix: MatrixPayload | None = None
    chronological_age: list[float] | None = None


class ResilienceCSDRequest(BaseModel):
    dataset_id: str | None = None
    data: MatrixPayload | None = None
    biomarker_cols: list[str] | None = None
    age_col: str = "age"
    n_strata: int = Field(6, ge=2, le=20)
    longitudinal: bool = False


class ResilienceRecoveryRequest(BaseModel):
    series: list[float] = Field(..., min_length=4)


class ErrorBody(BaseModel):
    code: str
    message: str
    detail: Any | None = None


class ErrorEnvelope(BaseModel):
    error: ErrorBody
