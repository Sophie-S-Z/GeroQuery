"""Wrappers around real published aging clocks from ``biolearn``.

We wrap, never rebuild. Horvath, Hannum, PhenoAge, GrimAge and friends have
published coefficients and reference implementations; re-deriving them here would
add a second, unvalidated copy of numbers that already exist.

Nothing in this module fabricates a clock. If ``biolearn`` is not installed,
:func:`library_clocks` returns an empty dict and the registry contains only the
transparent reference clocks.

Two things this layer is actually responsible for:

**Preserving ``predicted_outcome``.** biolearn's gallery contains far more than
age predictors — ``Smoking Status``, ``BMI``, ``Total Cholesterol``,
``Mortality Risk``. Flattening all of those to "aging clock" is the exact error
the PRD calls out: a clock's chronological-age accuracy says nothing about its
mortality prediction, and a smoking predictor is not an aging clock at all. Each
model's declared output is mapped explicitly, and anything unrecognized is passed
through under its own label rather than defaulted to ``chronological_age``.

**Orientation of the feature matrix.** GeroQuery passes samples as *rows*;
biolearn expects features as rows and samples as columns. The transpose happens
here so callers never have to know.

pyaging is deliberately *not* wrapped. It requires torch and ships its clocks as
downloaded model artifacts, so a wrapper that cannot be exercised in CI would be
untested code claiming real-clock support. Adding it is tracked work, not a
silent gap.
"""

from __future__ import annotations

import importlib.util
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from ..exceptions import ClockInputError
from ..models import ClockInfo

BIOLEARN = "biolearn"

# biolearn's `output` field -> GeroQuery's predicted_outcome vocabulary.
# Anything not listed is preserved verbatim (slugified) rather than coerced.
OUTCOME_MAP: dict[str, str] = {
    "Age (Years)": "chronological_age",
    "Age (years)": "chronological_age",
    "Age (days)": "chronological_age",
    "Age score": "chronological_age",
    "Human Cortex Age (Years)": "chronological_age",
    "Gestational Age": "gestational_age",
    "Gestational age": "gestational_age",
    "Mortality Risk": "mortality",
    "Mortality Risk by Organ": "mortality",
    "Mortality Adjusted Age (Years)": "mortality_adjusted_age",
    "Aging Rate (Years/Year)": "pace_of_aging",
    "Mitotic Age (Cell Divisions)": "mitotic_age",
    "Stem Cell Division Rate": "mitotic_age",
    "Telomere Length": "telomere_length",
}

# Units per outcome family, for the ones we recognize.
UNITS_MAP: dict[str, str] = {
    "chronological_age": "years",
    "gestational_age": "weeks",
    "mortality": "risk_score",
    "mortality_adjusted_age": "years",
    "pace_of_aging": "years_per_year",
    "mitotic_age": "cell_divisions",
    "telomere_length": "kb",
}

# biolearn model type -> the GeoData attribute it reads, and our input_type.
MODEL_TYPE_INPUT: dict[str, tuple[str, str]] = {
    "LinearMethylationModel": ("dnam", "methylation"),
    "GrimageModel": ("dnam", "methylation"),
    "SexEstimationModel": ("dnam", "methylation"),
    "DeconvolutionModel": ("dnam", "methylation"),
    "EpiTOC2Model": ("dnam", "methylation"),
    "AltumAgeModel": ("dnam", "methylation"),
    "PCLinearTransformationModel": ("dnam", "methylation"),
    "MiAgeModel": ("dnam", "methylation"),
    "GPAgeModel": ("dnam", "methylation"),
    "LinearTranscriptomicModel": ("rna", "expression"),
}

# Full CpG lists run to tens of thousands of entries. Exposing all of them in
# every /v1/clocks response would dwarf the rest of the payload.
MAX_ADVERTISED_FEATURES = 50


def _slug(text: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in text.strip().lower()).strip("_")


def map_outcome(output: str | None) -> str:
    """Map a biolearn ``output`` string to a predicted_outcome label.

    Unknown outputs are slugified and returned as-is. Defaulting them to
    ``chronological_age`` would label a BMI or smoking predictor an age clock.
    """
    if not output:
        return "unknown"
    if output in OUTCOME_MAP:
        return OUTCOME_MAP[output]
    return _slug(output)


@dataclass(frozen=True)
class LibraryClock:
    """Adapts one biolearn model to the registry's ``predict`` contract."""

    info: ClockInfo
    _model: Any
    _data_attr: str
    _all_features: tuple[str, ...] = ()

    @property
    def required_feature_set(self) -> tuple[str, ...]:
        """Every feature the model needs, not just the advertised subset."""
        return self._all_features or self.info.required_features

    def predict(self, matrix: pd.DataFrame, metadata: pd.DataFrame | None = None) -> np.ndarray:
        """Predict from a samples-as-rows matrix.

        Args:
            metadata: per-sample covariates indexed like ``matrix``. Some real
                clocks need them — GrimAge is trained on age and sex alongside
                methylation and raises without them. Passing an empty frame is
                what made GrimAge unusable in the first version of this wrapper.
        """
        required = set(self.required_feature_set)
        if required:
            missing = sorted(required - set(matrix.columns))
            if missing:
                raise ClockInputError(
                    f"Clock {self.info.clock_id!r} requires features not present in the input.",
                    detail={
                        "clock_id": self.info.clock_id,
                        "n_missing": len(missing),
                        "missing_features": missing[:MAX_ADVERTISED_FEATURES],
                        "n_required": len(required),
                        "provided": len(matrix.columns),
                    },
                )

        from biolearn.data_library import GeoData

        # biolearn wants features x samples; GeroQuery hands us samples x features.
        # Rebuilt from an owned ndarray rather than `.T`, for two reasons found by
        # running a real clock rather than a stub:
        #   1. biolearn's LinearModel.predict does `matrix_data.loc["intercept"] = 1`,
        #      so a view would mutate the caller's DataFrame.
        #   2. DunedinPACE quantile-normalizes in place and needs a writable
        #      buffer, else "assignment destination is read-only".
        # (2) only holds under pandas 2 — pandas 3 makes copy-on-write mandatory
        # and hands back read-only arrays regardless, which is why the `clocks`
        # extra pins `pandas<3`. Nothing outside biolearn can work around that.
        transposed = pd.DataFrame(
            matrix.to_numpy(dtype=float, copy=True).T,
            index=matrix.columns,
            columns=matrix.index,
        )
        meta = (
            pd.DataFrame(index=matrix.index) if metadata is None else metadata.reindex(matrix.index)
        )
        geo = GeoData(metadata=meta, **{self._data_attr: transposed})

        try:
            result = self._model.predict(geo)
        except ValueError as exc:
            raise ClockInputError(
                f"Clock {self.info.clock_id!r} could not run on this input: {exc}",
                detail={"clock_id": self.info.clock_id},
            ) from exc

        return _to_array(result, n=len(matrix))


def _to_array(result: Any, n: int) -> np.ndarray:
    """Normalize biolearn's return (DataFrame or Series) to a 1-D float array."""
    if isinstance(result, pd.DataFrame):
        # Single-output models return one column; take it rather than guessing a name.
        series = result.iloc[:, 0] if result.shape[1] >= 1 else pd.Series(dtype=float)
    elif isinstance(result, pd.Series):
        series = result
    else:
        series = pd.Series(result)
    values = np.asarray(series, dtype=float).ravel()
    if values.size != n:
        raise ClockInputError(
            "Clock returned a different number of predictions than samples provided.",
            detail={"expected": n, "got": int(values.size)},
        )
    return values


def _required_features(model: Any) -> tuple[str, ...]:
    """Best-effort feature list. Model types differ; absence is not an error."""
    getter = getattr(model, "methylation_sites", None)
    if callable(getter):
        try:
            return tuple(sorted(getter()))
        except Exception:
            return ()
    coefficients = getattr(model, "coefficients", None)
    if coefficients is not None and hasattr(coefficients, "index"):
        return tuple(sorted(str(i) for i in coefficients.index if str(i) != "intercept"))
    return ()


def _build_clock(name: str, definition: dict, gallery: Any) -> LibraryClock | None:
    model_type = (definition.get("model") or {}).get("type")
    mapping = MODEL_TYPE_INPUT.get(model_type) if isinstance(model_type, str) else None
    if mapping is None:
        # An unrecognized model type means we cannot know which GeoData slot it
        # reads. Skip it rather than register a clock that raises when called.
        return None
    data_attr, input_type = mapping

    try:
        model = gallery.get(name)
    except Exception:
        # A model whose coefficient file will not load is simply unavailable.
        return None

    outcome = map_outcome(definition.get("output"))
    features = _required_features(model)
    year = definition.get("year")
    tissue = definition.get("tissue") or "unspecified"
    species = definition.get("species") or "unspecified"

    return LibraryClock(
        info=ClockInfo(
            clock_id=f"biolearn:{name}",
            name=f"{name}{f' ({year})' if year else ''}",
            library=BIOLEARN,
            predicted_outcome=outcome,
            training_population=f"{species}, {tissue}" + (f", published {year}" if year else ""),
            input_type=input_type,
            units=UNITS_MAP.get(outcome, definition.get("output") or "unknown"),
            required_features=features[:MAX_ADVERTISED_FEATURES],
            notes=(
                f"Wrapped biolearn model. Declared output: {definition.get('output')!r}. "
                f"Requires {len(features)} features"
                + (
                    f" (first {MAX_ADVERTISED_FEATURES} listed)."
                    if len(features) > MAX_ADVERTISED_FEATURES
                    else "."
                )
                + (f" Source: {definition['source']}" if definition.get("source") else "")
            ),
        ),
        _model=model,
        _data_attr=data_attr,
        _all_features=features,
    )


def biolearn_available() -> bool:
    return importlib.util.find_spec(BIOLEARN) is not None


@dataclass(frozen=True)
class LibraryStatus:
    """Why the library-clock tier looks the way it does."""

    installed: bool
    usable: bool
    n_clocks: int
    n_skipped: int
    reason: str | None = None

    def to_dict(self) -> dict:
        return {
            "installed": self.installed,
            "usable": self.usable,
            "n_clocks": self.n_clocks,
            "n_skipped": self.n_skipped,
            "reason": self.reason,
        }


def load_library_clocks() -> tuple[dict[str, LibraryClock], LibraryStatus]:
    """Wrappable biolearn models plus a diagnosis of anything that did not load.

    The status is returned, not just logged, because "biolearn is installed but
    the registry has zero real clocks" is indistinguishable from "biolearn is not
    installed" at the ``/v1/clocks`` boundary — and an earlier version of this
    function swallowed the exception, so a missing transitive dependency
    (biolearn imports ``torch`` at module scope without declaring it) presented
    as a silently empty clock list.
    """
    if not biolearn_available():
        return {}, LibraryStatus(
            installed=False,
            usable=False,
            n_clocks=0,
            n_skipped=0,
            reason="biolearn is not installed; install with `pip install -e '.[clocks]'`.",
        )
    try:
        from biolearn.model_gallery import ModelGallery

        gallery = ModelGallery()
        definitions = gallery.model_definitions
    except Exception as exc:
        return {}, LibraryStatus(
            installed=True,
            usable=False,
            n_clocks=0,
            n_skipped=0,
            reason=f"biolearn is installed but its gallery could not load: {exc!r}",
        )

    out: dict[str, LibraryClock] = {}
    skipped: list[str] = []
    for name, definition in definitions.items():
        try:
            clock = _build_clock(name, definition, gallery)
        except Exception:
            clock = None
        if clock is None:
            skipped.append(name)
        else:
            out[clock.info.clock_id] = clock

    reason = None
    if not out:
        reason = (
            f"biolearn loaded but none of its {len(definitions)} models could be "
            "wrapped (unknown model types, or coefficient files unavailable)."
        )
    elif skipped:
        reason = f"{len(skipped)} model(s) skipped: {', '.join(sorted(skipped)[:8])}"

    return out, LibraryStatus(
        installed=True,
        usable=bool(out),
        n_clocks=len(out),
        n_skipped=len(skipped),
        reason=reason,
    )


def library_clocks() -> dict[str, LibraryClock]:
    """Wrappable biolearn models. See :func:`load_library_clocks` for the reason
    an empty result is empty."""
    return load_library_clocks()[0]
