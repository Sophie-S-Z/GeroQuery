"""Perturb-seq — CRISPRi knockdown screens read out across the transcriptome.

Replogle et al. (*Cell* 2022) knocked down one gene at a time in human cell
lines and sequenced the result, giving, for each knockdown, the shift in every
measured gene. That is the one public dataset that answers "what happens
downstream when this gene is switched off" at scale, and it is CC BY 4.0.

**What is read.** Only the pseudobulk files: the per-knockdown average over
cells. Every number this layer computes is such an average, so the 8.7 GB and
66 GB single-cell versions of the same screens would cost three orders of
magnitude more to reach the same answer.

**Two screens, deliberately.** RPE1 (retinal pigment epithelium, not a cancer
line) and K562 (chronic myeloid leukaemia, screened genome-wide). Neither is
aged tissue and both are immortalized dividing cells, so no single screen can
support a claim about aging. Agreement between two unrelated cell lines is the
only replication available, and the exporter requires it — the same standard
the panel already applies to microarray against RNA-seq.

**Selection is a rule, not a list.** A knockdown is scored when:

* its response is finite in every readout gene kept (two readout genes carry
  infinities in the RPE1 release and are dropped for every perturbation, so the
  gene set stays identical across knockdowns);
* the screen's own energy test rejects "this knockdown changed nothing" at
  :data:`EFFECT_P`. Replogle computed that test per perturbation; 2,289 of 2,679
  RPE1 knockdowns and 5,401 of 11,258 K562 knockdowns pass. Scoring the rest
  would be correlating noise with the aging direction and reporting the result.

**The built-in positive control.** A knockdown must lower its own target. Across
the screens it does, at a median z of -0.64 against an off-target spread of
0.29, in 98.8% of perturbations. :func:`on_target` recomputes it, and the test
suite pins it: when the headline results of this layer are negative, the control
is what separates "measured nothing" from "loaded the file wrong".
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ..exceptions import SourceError

if TYPE_CHECKING:
    import h5py

# The screen's own test that a knockdown did anything at all.
EFFECT_P = 0.05

# Rows whose guides target no gene: the screen's own negative control.
CONTROL_SYMBOL = "non-targeting"

# Manifest keys, with the cell line each was run in.
SCREENS: dict[str, str] = {
    "replogle_rpe1": "RPE1",
    "replogle_k562_gwps": "K562",
}


@dataclass(frozen=True)
class Screen:
    """One pseudobulk screen, after the selection rule."""

    key: str
    cell_line: str
    responses: np.ndarray
    """(perturbations x readout genes) normalized expression shift."""

    readout_genes: np.ndarray
    """Ensembl id per column."""

    readout_symbols: np.ndarray
    """Symbol per column, for display only; joins go through the Ensembl id."""

    perturbed_genes: np.ndarray
    """Ensembl id of the gene knocked down, per row."""

    symbols: np.ndarray
    """Symbol of the gene knocked down, per row."""

    has_effect: np.ndarray
    """Per row: did the screen's own energy test find any effect?"""

    is_control: np.ndarray
    """Per row: is this a non-targeting control rather than a knockdown?"""

    def __post_init__(self) -> None:
        rows, cols = self.responses.shape
        for name, got, want in (
            ("readout_genes", len(self.readout_genes), cols),
            ("readout_symbols", len(self.readout_symbols), cols),
            ("perturbed_genes", len(self.perturbed_genes), rows),
            ("symbols", len(self.symbols), rows),
            ("has_effect", len(self.has_effect), rows),
            ("is_control", len(self.is_control), rows),
        ):
            if got != want:
                raise SourceError(f"{self.key}: {name} has {got} entries, expected {want}.")

    @property
    def scorable(self) -> np.ndarray:
        """Rows that are a real knockdown, had an effect, and can be joined.

        Controls target nothing, unmapped targets cannot be joined to the panel,
        and a knockdown the energy test cannot distinguish from no change would
        contribute noise dressed as a measurement.
        """
        return self.has_effect & ~self.is_control & (self.perturbed_genes != "")

    @property
    def n_scored(self) -> int:
        return int(self.scorable.sum())


def _text(values: np.ndarray) -> np.ndarray:
    """Decode an HDF5 byte-string column."""
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in values])


def _gene_symbols(handle: h5py.File) -> np.ndarray:
    """Readout symbols, stored as codes into a category table."""
    codes = handle["var/gene_name"][:]
    if "var/__categories/gene_name" not in handle:
        return _text(codes)
    categories = _text(handle["var/__categories/gene_name"][:])
    return categories[codes]


def read_screen(path: Path | str, key: str) -> Screen:
    """Load one pinned pseudobulk file and apply the selection rule.

    The perturbation index is ``<id>_<SYMBOL>_<protospacers>_<ENSG>``, which
    carries both the symbol and the Ensembl id of the gene knocked down. The
    Ensembl id is what the rest of GeroQuery joins on, so it is the identity
    kept; the symbol is carried for display only.
    """
    cell_line = SCREENS.get(key)
    if cell_line is None:
        raise SourceError(f"Unknown Perturb-seq screen {key!r}. Known: {sorted(SCREENS)}.")
    # Imported here, as in archs4: only reading a screen file needs HDF5, so the
    # API and the exporter's scoring can run without the optional package.
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise SourceError(
            "Reading a Perturb-seq screen needs the optional 'h5py' package: "
            "pip install 'geroquery[perturbation]'."
        ) from exc
    with h5py.File(path, "r") as handle:
        responses = handle["X"][:]
        readout = _text(handle["var/gene_id"][:])
        readout_symbols = _gene_symbols(handle)
        index = _text(handle["obs/gene_transcript"][:])
        effect_p = handle["obs/energy_test_p_value"][:]

    finite = np.isfinite(responses).all(axis=0)
    if not finite.any():
        raise SourceError(f"{key}: no readout gene is finite across every perturbation.")
    responses, readout = responses[:, finite], readout[finite]
    readout_symbols = readout_symbols[finite]

    parts = [name.split("_") for name in index]
    if any(len(p) < 4 for p in parts):
        raise SourceError(f"{key}: perturbation names are not '<id>_<symbol>_<guides>_<gene>'.")
    symbols = np.array([p[1] for p in parts])
    is_control = symbols == CONTROL_SYMBOL
    # Unmapped targets carry the literal "nan" and cannot be joined to anything
    # GeroQuery knows; controls target no gene at all. Both keep their row and
    # lose their identity, so neither can be mistaken for a scored knockdown.
    perturbed = np.array(
        [p[-1] if p[-1].startswith("ENSG") else "" for p in parts],
        dtype=object,
    )
    if not (perturbed != "").any():
        raise SourceError(f"{key}: no perturbation carries an Ensembl gene id.")

    return Screen(
        key=key,
        cell_line=cell_line,
        responses=responses,
        readout_genes=readout,
        readout_symbols=readout_symbols,
        perturbed_genes=perturbed.astype(str),
        symbols=symbols,
        has_effect=np.asarray(effect_p) < EFFECT_P,
        is_control=is_control,
    )


def symbol_map(screen: Screen) -> dict[str, str]:
    """Ensembl id to symbol for every gene the screen touches, read out or hit.

    Display only. The panel joins on the Ensembl id because the screens' symbols
    and the panel's drift apart, but a table of bare accessions is unreadable.
    """
    mapping = dict(zip(screen.readout_genes, screen.readout_symbols, strict=True))
    for gene, symbol, control in zip(
        screen.perturbed_genes, screen.symbols, screen.is_control, strict=True
    ):
        if gene and not control:
            mapping[gene] = symbol
    return mapping


def on_target(screen: Screen) -> np.ndarray:
    """Each knockdown's effect on its own target, where that gene is read out.

    The positive control. Returns the response of the perturbed gene itself, one
    value per perturbation that has one, which should be negative: knocking a
    gene down lowers it.
    """
    column = {gene: i for i, gene in enumerate(screen.readout_genes)}
    return np.array(
        [
            screen.responses[row, column[gene]]
            for row, gene in enumerate(screen.perturbed_genes)
            if screen.scorable[row] and gene in column
        ]
    )
