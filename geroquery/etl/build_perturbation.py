"""Build the perturbation layer: knockdown screens against the aging direction.

Run: ``python -m geroquery.etl.build_perturbation``

Answers three questions from the pinned Perturb-seq screens, and is explicit
about which of them the data can carry:

1. **Which aging genes are testable at all.** Of the genes GeroQuery resolves a
   direction for, which have been knocked down in a human cell line, and which
   are only read out, or absent. This is a coverage fact and needs no model.
2. **Does knocking a gene down move the cell toward or away from the aging
   signature.** Scored by :mod:`geroquery.harmonize.perturbation`, which
   removes the damage confound. A ranking, not a verdict: no gene reaches an
   FDR-significant alignment in either screen, so the exporter ships the score
   with its cross-screen agreement and flags only the genes whose extremes
   replicate.
3. **Whether pooled aging claims hold under perturbation.** A corpus-level
   association between a gene's own aging direction and its knockdown's
   alignment. Measured and reported whichever way it comes out; it comes out
   near zero.

Output lands in ``frontend/public/data/``:

``perturbation.parquet``
    One row per human gene the screens touch, with each screen's alignment,
    adjusted score and coverage flags, plus the cross-screen verdict.

The headline numbers, the positive control, and the confound diagnostics go
into ``meta.json`` under ``perturbation`` so the interface and the test suite
read the same figures.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .. import config
from ..harmonize.perturbation import adjust_for_magnitude, align_many
from ..sources.fetch import fetch_artifact
from ..sources.manifest import PERTURB_SEQ
from ..sources.perturbseq import SCREENS, Screen, on_target, read_screen, symbol_map

DEFAULT_OUT = config.FRONTEND_DATA

# A gene's knockdown replicates when both screens rank it this extreme and
# agree on the sign. The full ranking agrees at 52.7%, so only the tails are
# worth reporting; the top decile is where the cross-screen overlap runs 2-4x
# above chance.
EXTREME_QUANTILE = 0.10

# Below this many shared readout genes the alignment is over too thin a slice
# of the aging direction to mean anything.
MIN_AGING_GENES = 200


@dataclass(frozen=True)
class ScoredScreen:
    """One screen's scores, keyed by the Ensembl id of the perturbed gene."""

    key: str
    cell_line: str
    scores: pd.DataFrame
    readout_genes: set[str]
    attempted_genes: set[str]
    """Genes knocked down at all, whether or not the knockdown did anything."""

    symbols: dict[str, str]
    """Ensembl id to symbol for every gene this screen touches."""

    n_perturbations: int
    on_target_median: float
    on_target_fraction_down: float
    magnitude_confound: float
    n_aging_genes: int
    n_controls: int
    control_alignment_sd: float
    knockdown_alignment_sd: float


def aging_direction(pooled: pd.DataFrame) -> pd.Series:
    """Pooled human effect per gene, over genes GeroQuery resolves a direction for.

    Restricting to resolved genes keeps the direction made of estimates the
    panel is willing to state. Using the shrunk estimate instead changes the
    resulting scores by r = 0.9995, so the published ``g`` is used.
    """
    human = pooled[pooled["species"] == "human"]
    resolved = human[human["verdict"] != "no_evidence"]
    return pd.Series(
        resolved["g"].to_numpy(dtype=float), index=resolved["gene_id"].astype(str), name="g"
    )


def score_screen(screen: Screen, direction: pd.Series) -> ScoredScreen:
    """Score every knockdown in one screen against the aging direction."""
    shared = [g for g in screen.readout_genes if g in direction.index]
    if len(shared) < MIN_AGING_GENES:
        raise ValueError(
            f"{screen.key}: only {len(shared)} readout genes carry an aging direction, "
            f"need {MIN_AGING_GENES}."
        )
    column = {gene: i for i, gene in enumerate(screen.readout_genes)}
    take = np.array([column[g] for g in shared])
    responses = screen.responses[:, take]
    position = {gene: i for i, gene in enumerate(shared)}

    rows = np.flatnonzero(screen.scorable)
    skips = [position.get(screen.perturbed_genes[r]) for r in rows]
    scored = align_many(responses[rows], direction.loc[shared].to_numpy(), skip=skips)

    frame = pd.DataFrame(
        {
            "gene_id": screen.perturbed_genes[rows],
            "symbol": screen.symbols[rows],
            "alignment": [s.alignment for s in scored],
            "magnitude": [s.magnitude for s in scored],
        }
    )
    # A flat response has no direction, and a zero magnitude has no logarithm.
    frame = frame[np.isfinite(frame["alignment"]) & (frame["magnitude"] > 0)]
    frame["z"] = adjust_for_magnitude(frame["alignment"], frame["magnitude"])

    attempted = {
        g for g, c in zip(screen.perturbed_genes, screen.is_control, strict=True) if g and not c
    }
    control = on_target(screen)
    # The screen's own negative control, scored the same way: guides that target
    # no gene. Their alignment spread is what this statistic reads off a cell
    # nothing was done to, and real knockdowns must spread wider to mean anything.
    control_rows = np.flatnonzero(screen.is_control)
    control_scores = align_many(
        responses[control_rows], direction.loc[shared].to_numpy(), skip=None
    )
    control_alignment = np.array([s.alignment for s in control_scores if np.isfinite(s.alignment)])
    # Several guide sets can target one gene; a gene gets one score, their mean.
    per_gene = (
        frame.groupby("gene_id")
        .agg(symbol=("symbol", "first"), alignment=("alignment", "mean"), z=("z", "mean"))
        .reset_index()
    )
    return ScoredScreen(
        key=screen.key,
        cell_line=screen.cell_line,
        scores=per_gene.set_index("gene_id"),
        readout_genes=set(screen.readout_genes),
        attempted_genes=attempted,
        symbols=symbol_map(screen),
        n_perturbations=len(screen.perturbed_genes),
        on_target_median=float(np.median(control)),
        on_target_fraction_down=float((control < 0).mean()),
        magnitude_confound=float(np.corrcoef(frame["alignment"], np.log(frame["magnitude"]))[0, 1]),
        n_aging_genes=len(shared),
        n_controls=len(control_alignment),
        control_alignment_sd=float(control_alignment.std()) if len(control_alignment) else 0.0,
        knockdown_alignment_sd=float(frame["alignment"].std()),
    )


def combine(scored: list[ScoredScreen]) -> pd.DataFrame:
    """One row per gene the screens touch, with the cross-screen verdict.

    ``replicates`` requires both screens to rank the gene in the extreme
    :data:`EXTREME_QUANTILE` tail *and* to agree on the sign. Everything else is
    reported as scored-in-one, discordant, or read out but never perturbed.
    """
    genes: set[str] = set()
    for s in scored:
        genes |= set(s.scores.index) | s.readout_genes | s.attempted_genes
    frame = pd.DataFrame({"gene_id": sorted(genes)}).set_index("gene_id")

    # Symbols come from every gene a screen touches, not only the scored ones,
    # or four rows in five would show a bare accession.
    names: dict[str, str] = {}
    for s in reversed(scored):
        names.update(s.symbols)
    frame["symbol"] = frame.index.map(lambda gene: names.get(gene, ""))
    for s in scored:
        prefix = s.cell_line.lower()
        joined = s.scores.reindex(frame.index)
        frame[f"{prefix}_alignment"] = joined["alignment"].round(4)
        frame[f"{prefix}_z"] = joined["z"].round(4)
        frame[f"{prefix}_scored"] = joined["z"].notna()
        frame[f"{prefix}_readout"] = frame.index.isin(s.readout_genes)
        frame[f"{prefix}_perturbed"] = frame.index.isin(s.attempted_genes)

    lines = [s.cell_line.lower() for s in scored]
    both = frame[[f"{c}_scored" for c in lines]].all(axis=1)
    z = frame[[f"{c}_z" for c in lines]]
    agree = both & (np.sign(z[f"{lines[0]}_z"]) == np.sign(z[f"{lines[1]}_z"]))

    extreme = pd.Series(True, index=frame.index)
    for s, c in zip(scored, lines, strict=True):
        cut = s.scores["z"].abs().quantile(1 - EXTREME_QUANTILE)
        extreme &= frame[f"{c}_z"].abs() >= cut
    frame["replicates"] = (agree & extreme).fillna(False)

    direction = np.where(z[f"{lines[0]}_z"] > 0, "toward_aging", "away_from_aging")
    frame["direction"] = np.where(
        frame["replicates"], direction, np.where(both & ~agree, "discordant", "")
    )
    frame["n_screens"] = frame[[f"{c}_scored" for c in lines]].sum(axis=1).astype(int)
    return frame.reset_index()


def _coverage(combined: pd.DataFrame, pooled: pd.DataFrame, lines: list[str]) -> dict:
    """What the screens can and cannot be asked about GeroQuery's aging genes."""
    human = pooled[pooled["species"] == "human"]
    resolved = set(human.loc[human["verdict"] != "no_evidence", "gene_id"].astype(str))
    scored = set(combined.loc[combined["n_screens"] > 0, "gene_id"])
    attempted = set(
        combined.loc[combined[[f"{c}_perturbed" for c in lines]].any(axis=1), "gene_id"]
    )
    readout = set(combined.loc[combined[[f"{c}_readout" for c in lines]].any(axis=1), "gene_id"])
    return {
        "resolved_genes": len(resolved),
        # Knocked down, and the screen's energy test saw the cell respond.
        "resolved_scored": len(resolved & scored),
        # Knocked down, and nothing measurable happened: tested, not untested.
        "resolved_no_effect": len((resolved & attempted) - scored),
        # Measured as a readout gene, but never itself knocked down.
        "resolved_readout_only": len((resolved & readout) - attempted),
        "resolved_absent": len(resolved - attempted - readout),
    }


def _cross_modality(combined: pd.DataFrame, pooled: pd.DataFrame, lines: list[str]) -> dict:
    """Does a gene's own aging direction predict its knockdown's alignment?

    If genes that rise with age were driving the aging transcriptome, knocking
    them down should move cells away from it, giving a negative association.
    Reported as measured.
    """
    human = pooled[pooled["species"] == "human"].set_index("gene_id")
    out = {}
    for line in lines:
        scored = combined[combined[f"{line}_scored"]].set_index("gene_id")
        shared = scored.index.intersection(human.index[human["verdict"] != "no_evidence"])
        if len(shared) < MIN_AGING_GENES:  # pragma: no cover - guarded upstream
            continue
        rho, p = stats.spearmanr(
            human.loc[shared, "g"].to_numpy(dtype=float), scored.loc[shared, f"{line}_z"]
        )
        out[line] = {"n": len(shared), "spearman": round(float(rho), 4), "p": float(p)}
    return out


def build(out_dir: Path | None = None, *, allow_network: bool | None = None) -> dict:
    """Score both screens, write the artifact, and return the meta block."""
    out_dir = out_dir or DEFAULT_OUT
    pooled = pd.read_parquet(out_dir / "pooled.parquet")
    direction = aging_direction(pooled)

    scored = []
    for key in SCREENS:
        path = fetch_artifact(key, allow_network=allow_network)
        scored.append(score_screen(read_screen(path, key), direction))

    combined = combine(scored)
    combined.to_parquet(out_dir / "perturbation.parquet", index=False, compression="snappy")

    lines = [s.cell_line.lower() for s in scored]
    both = combined[combined["n_screens"] == len(scored)]
    z0, z1 = both[f"{lines[0]}_z"], both[f"{lines[1]}_z"]
    meta = {
        "release": "Replogle et al. 2022 processed Perturb-seq (pseudobulk)",
        "license": PERTURB_SEQ["replogle_rpe1"].license,
        "attribution": PERTURB_SEQ["replogle_rpe1"].attribution,
        "screens": [
            {
                "key": s.key,
                "cell_line": s.cell_line,
                "n_perturbations": s.n_perturbations,
                "n_scored": int(len(s.scores)),
                "n_aging_genes": s.n_aging_genes,
                "on_target_median": round(s.on_target_median, 4),
                "on_target_fraction_down": round(s.on_target_fraction_down, 4),
                "magnitude_confound": round(s.magnitude_confound, 4),
                "n_controls": s.n_controls,
                "control_alignment_sd": round(s.control_alignment_sd, 4),
                "knockdown_alignment_sd": round(s.knockdown_alignment_sd, 4),
            }
            for s in scored
        ],
        "n_scored_in_both": int(len(both)),
        "cross_screen_spearman": round(float(stats.spearmanr(z0, z1).statistic), 4),
        "cross_screen_sign_agreement": round(float((np.sign(z0) == np.sign(z1)).mean()), 4),
        "n_replicating": int(combined["replicates"].sum()),
        "extreme_quantile": EXTREME_QUANTILE,
        "coverage": _coverage(combined, pooled, lines),
        "cross_modality": _cross_modality(combined, pooled, lines),
    }

    meta_path = out_dir / "meta.json"
    payload = json.loads(meta_path.read_text(encoding="utf-8"))
    payload["perturbation"] = meta
    meta_path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    return meta


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args(argv)
    meta = build(args.out_dir)
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
