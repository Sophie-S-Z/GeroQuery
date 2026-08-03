"""Wrappers around real published aging clocks from ``pyaging``.

Second library tier, alongside :mod:`geroquery.clocks.library` (biolearn). Between
them the registry exposes 200+ real published clocks without reimplementing any.

pyaging works differently from biolearn in one way that shapes this whole module:

**Artifacts are downloaded, not shipped.** biolearn ships coefficient CSVs inside
its wheel, so its clocks can be enumerated *and* have their feature lists read
offline. pyaging keeps a catalog and one model artifact per clock on Hugging Face,
fetched on demand. So:

* **Registration is metadata-only.** We list clocks from the cached catalog and do
  not know their required features until the artifact is downloaded. Eagerly
  downloading 173 artifacts to populate `required_features` would mean a
  multi-hundred-megabyte fetch just to answer ``GET /v1/clocks``.
* **``required_features`` is therefore empty for pyaging clocks**, and the
  wrapper's own pre-flight check is skipped — pyaging reports missing features
  itself, and we translate that into a :class:`ClockInputError`.
* **The first prediction for a given clock hits the network.** Subsequent ones use
  the Hugging Face cache.

Registration is skipped entirely when the catalog is not cached and network access
is off, so an offline run is never blocked by this tier.
"""

from __future__ import annotations

import importlib.util
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .. import config
from ..exceptions import ClockInputError
from ..models import ClockInfo
from .library import LibraryStatus

PYAGING = "pyaging"

# pyaging's `predicts` entries -> GeroQuery's predicted_outcome vocabulary.
# Anything unlisted is slugified and kept, never defaulted to an age outcome.
OUTCOME_MAP: dict[str, str] = {
    "chronological age": "chronological_age",
    "transcriptomic age": "chronological_age",
    "biological age": "biological_age",
    "phenotypic age": "biological_age",
    "physical-fitness biological age": "biological_age",
    "gestational age": "gestational_age",
    "mortality risk": "mortality",
    "pace of aging": "pace_of_aging",
    "mitotic age": "mitotic_age",
    "replicative history": "mitotic_age",
    "cell-population passage age": "mitotic_age",
    "leukocyte telomere length": "telomere_length",
    "intervention-responsive epigenetic age": "biological_age",
}

UNITS_MAP: dict[str, str] = {
    "chronological_age": "years",
    "biological_age": "years",
    "gestational_age": "weeks",
    "mortality": "risk_score",
    "pace_of_aging": "years_per_year",
    "mitotic_age": "cell_divisions",
    "telomere_length": "kb",
}

# pyaging silently imputes absent features from reference values. Past this
# fraction the prediction says more about pyaging's reference cohort than about
# the samples, so we refuse rather than return a confident-looking number.
MAX_IMPUTED_FRACTION = 0.5

# pyaging `data_type` -> GeroQuery input_type.
DATA_TYPE_INPUT: dict[str, str] = {
    "DNA methylation": "methylation",
    "transcriptomics": "expression",
    "clinical biomarkers": "clinical",
    "histone modification": "histone_modification",
    "chromatin accessibility": "chromatin_accessibility",
}


def _slug(text: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in str(text).strip().lower()).strip("_")


def map_outcome(predicts: Any) -> str:
    """Map pyaging's ``predicts`` list to one predicted_outcome label.

    A clock predicting "smoking exposure" or "grip strength" must not come back
    as an age clock, so unmapped values are slugified and preserved.
    """
    if isinstance(predicts, (list, tuple)):
        first = predicts[0] if predicts else None
    else:
        first = predicts
    if not first:
        return "unknown"
    key = str(first).strip().lower()
    return OUTCOME_MAP.get(key, _slug(key))


def pyaging_available() -> bool:
    return importlib.util.find_spec(PYAGING) is not None


@dataclass(frozen=True)
class PyagingClock:
    """Adapts one pyaging clock to the registry's ``predict`` contract."""

    info: ClockInfo
    clock_name: str

    @property
    def required_feature_set(self) -> tuple[str, ...]:
        """Empty by design — see the module docstring. pyaging validates its own
        feature requirements once the artifact is downloaded."""
        return ()

    def predict(self, matrix: pd.DataFrame, metadata: pd.DataFrame | None = None) -> np.ndarray:
        """Predict from a samples-as-rows matrix.

        pyaging takes samples as rows, matching GeroQuery's orientation, so unlike
        the biolearn wrapper no transpose is needed. The frame is still copied:
        ``df_to_adata`` and the imputation step both write into what they are given.
        """
        import pyaging as pya

        frame = matrix.copy()
        adata = None
        try:
            adata = pya.preprocess.df_to_adata(frame, verbose=False)
            pya.pred.predict_age(adata, self.clock_name, verbose=False)
        except NameError as exc:
            # pyaging raises a bare, message-less NameError when 100% of a
            # clock's features are absent. Propagating that verbatim tells the
            # caller nothing; it reads as an internal bug rather than "your
            # matrix has none of the CpGs this clock needs".
            raise ClockInputError(
                f"Clock {self.info.clock_id!r} found none of its required features in the "
                f"input. Check that the matrix columns are {self.info.input_type} feature "
                f"identifiers (e.g. CpG ids like 'cg00075967') rather than sample ids.",
                detail={
                    "clock_id": self.info.clock_id,
                    "clock_name": self.clock_name,
                    "provided_features": len(matrix.columns),
                    "example_provided": list(matrix.columns[:5]),
                },
            ) from exc
        except Exception as exc:
            raise ClockInputError(
                f"Clock {self.info.clock_id!r} could not run on this input: {exc}",
                detail={"clock_id": self.info.clock_id, "clock_name": self.clock_name},
            ) from exc

        # pyaging imputes absent features from reference values rather than
        # failing, so a prediction can be built from very little real data.
        # Surfacing the coverage is the difference between a usable number and a
        # confidently-reported imputation artifact.
        percent_na = adata.uns.get(f"{self.clock_name}_percent_na")
        if percent_na is not None and percent_na > MAX_IMPUTED_FRACTION * 100:
            raise ClockInputError(
                f"Clock {self.info.clock_id!r} is missing {percent_na:.1f}% of its features; "
                f"pyaging would impute them from reference values, making the prediction an "
                f"artifact of the reference rather than a measurement of these samples.",
                detail={
                    "clock_id": self.info.clock_id,
                    "percent_missing": float(percent_na),
                    "max_allowed_percent": MAX_IMPUTED_FRACTION * 100,
                },
            )

        column = self.clock_name
        if column not in adata.obs.columns:
            # pyaging names the output column after the clock; if that ever
            # changes, fail loudly rather than return the wrong column.
            raise ClockInputError(
                f"Clock {self.info.clock_id!r} produced no {column!r} column.",
                detail={"available": list(adata.obs.columns)},
            )
        values = np.asarray(adata.obs[column], dtype=float).ravel()
        if values.size != len(matrix):
            raise ClockInputError(
                "Clock returned a different number of predictions than samples provided.",
                detail={"expected": len(matrix), "got": int(values.size)},
            )
        return values


def _catalog(allow_network: bool) -> dict[str, dict]:
    """Load pyaging's clock catalog, preferring the local Hugging Face cache."""
    import os

    import pyaging as pya
    from pyaging.logger import LoggerManager

    previous = os.environ.get("HF_HUB_OFFLINE")
    if not allow_network:
        os.environ["HF_HUB_OFFLINE"] = "1"
    try:
        logger = LoggerManager.gen_logger("geroquery")
        # Exactly two positional args, and no more. pyaging's @progress decorator
        # does `logger = args[-1]`, so an all-keyword call raises IndexError and
        # passing indent_level positionally makes it treat the int as the logger.
        return pya.utils.load_clock_metadata("pyaging_data", logger)
    finally:
        if previous is None:
            os.environ.pop("HF_HUB_OFFLINE", None)
        else:
            os.environ["HF_HUB_OFFLINE"] = previous


def _build_clock(name: str, meta: dict) -> PyagingClock | None:
    data_type = meta.get("data_type")
    input_type = DATA_TYPE_INPUT.get(str(data_type), _slug(data_type or "unknown"))
    outcome = map_outcome(meta.get("predicts"))
    species = meta.get("species") or "unspecified"
    tissue = meta.get("tissue") or []
    if isinstance(tissue, (list, tuple)):
        tissue = ", ".join(str(t) for t in tissue) or "unspecified"
    year = meta.get("year")

    citation = meta.get("citation") or ""
    doi = meta.get("doi") or ""
    notes = meta.get("notes") or ""

    return PyagingClock(
        info=ClockInfo(
            clock_id=f"pyaging:{name}",
            name=f"{name}{f' ({year})' if year else ''}",
            library=PYAGING,
            predicted_outcome=outcome,
            training_population=f"{species}, {tissue}" + (f", published {year}" if year else ""),
            input_type=input_type,
            units=UNITS_MAP.get(outcome, "unknown"),
            # Deliberately empty: pyaging downloads the artifact lazily, so the
            # feature list is not knowable at registration time.
            required_features=(),
            notes=(
                f"Wrapped pyaging clock. Predicts: {meta.get('predicts')}. "
                f"Feature list resolved on first use (artifact downloaded from "
                f"Hugging Face). {notes} {citation} {doi}".strip()
            ),
        ),
        clock_name=name,
    )


def load_pyaging_clocks(
    allow_network: bool | None = None,
) -> tuple[dict[str, PyagingClock], LibraryStatus]:
    """Registered pyaging clocks plus a diagnosis of anything that did not load."""
    if not pyaging_available():
        return {}, LibraryStatus(
            installed=False,
            usable=False,
            n_clocks=0,
            n_skipped=0,
            reason="pyaging is not installed; install with `pip install -e '.[clocks]'`.",
        )

    allow_network = config.ALLOW_NETWORK if allow_network is None else allow_network
    try:
        catalog = _catalog(allow_network)
    except Exception as exc:
        return {}, LibraryStatus(
            installed=True,
            usable=False,
            n_clocks=0,
            n_skipped=0,
            reason=(
                "pyaging is installed but its clock catalog is unavailable "
                f"({'network disabled and nothing cached' if not allow_network else exc!r}). "
                "Run once with GEROQUERY_ALLOW_NETWORK=1 to populate the cache."
            ),
        )

    out: dict[str, PyagingClock] = {}
    skipped: list[str] = []
    for name, meta in catalog.items():
        try:
            clock = _build_clock(name, meta if isinstance(meta, dict) else {})
        except Exception:
            clock = None
        if clock is None:
            skipped.append(name)
        else:
            out[clock.info.clock_id] = clock

    return out, LibraryStatus(
        installed=True,
        usable=bool(out),
        n_clocks=len(out),
        n_skipped=len(skipped),
        reason=(f"{len(skipped)} clock(s) skipped" if skipped else None),
    )


def pyaging_clocks() -> dict[str, PyagingClock]:
    return load_pyaging_clocks()[0]
