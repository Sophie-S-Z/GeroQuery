"""What a knockdown predictor is worth once you measure its baseline.

Run: ``python -m geroquery.etl.perturbation_benchmark``

The "virtual cell" demo is: train on genes that were knocked down, hold one
out, predict what its knockdown does, and count how many of the genes that
moved most went the right way. The number comes out high. This module asks the
question that decides whether the number means anything — **how high is it for
a model that is not told which gene was knocked down?**

Most knockdowns of a consequential gene produce a similar response: the cell
stalls and stresses. A model that always predicts that one average profile
therefore gets most directions right on any held-out gene, having learned
nothing gene-specific. That is the control, and it is the one the demos omit.

Both models are evaluated under the same k-fold split, and the baseline is
built from the training folds only, so neither sees the held-out gene. The
ridge penalty is chosen by an inner split of the training folds, never on the
test fold: picking it per gene against the held-out answer adds about four
points of accuracy to a model with no gene-specific skill at all, which is
large enough to invent the very result this module exists to check.

This answers the fourth question in ``docs/PERTURBATION.md`` and deliberately
ships nothing into the product. A predictor that cannot beat "ignore the input"
has no business serving answers next to measured ones.

Features: a gene's own column in the screen — how it responds when *other*
genes are knocked down. That exists for every measured gene, including ones
never perturbed, which is what makes prediction for an unperturbed gene
possible at all. Target: that gene's row, what its own knockdown does.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy import stats

from ..sources.fetch import fetch_artifact
from ..sources.perturbseq import SCREENS, Screen, read_screen

# How many of the most-moved genes the direction is scored on, as the demos do.
TOP_N = 25

# Ridge penalties tried, as fractions of the largest squared singular value of
# the training features. Relative rather than absolute so the grid means the
# same thing on a screen of 5,000 genes and on a planted test of 160: an
# absolute grid tuned for one silently regularizes the other to a constant.
# The best per fold is reported, which favours the model.
RELATIVE_PENALTIES = (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0)

FOLDS = 5


@dataclass(frozen=True)
class Benchmark:
    """One screen's comparison. Accuracies are fractions of ``TOP_N``."""

    screen: str
    cell_line: str
    n_genes: int
    baseline_accuracy: float
    model_accuracy: float
    baseline_correlation: float
    model_correlation: float
    gain_pp: float
    gain_p: float

    @property
    def model_wins(self) -> bool:
        return self.gain_pp > 0.0


def _top_direction(predicted: np.ndarray, actual: np.ndarray) -> float:
    top = np.argsort(-np.abs(actual))[:TOP_N]
    return float((np.sign(predicted[top]) == np.sign(actual[top])).mean())


def _ridge_predictions(
    train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray
) -> list[np.ndarray]:
    """Ridge predictions for every penalty, from one SVD.

    Solving ``(X'X + lam I)^-1 X'Y`` per penalty would repeat the expensive
    part; the SVD makes each penalty a rescaling of the same factors and is
    exact, not an approximation.
    """
    u, s, vt = np.linalg.svd(train_x, full_matrices=False)
    uy = u.T @ train_y
    projected = test_x @ vt.T
    scale = float(s[0] ** 2)
    return [
        projected @ ((s / (s**2 + alpha * scale))[:, None] * uy) for alpha in RELATIVE_PENALTIES
    ]


def _choose_penalty(features: np.ndarray, targets: np.ndarray) -> int:
    """Index of the best penalty, measured on a split of the training genes only."""
    cut = max(1, len(targets) // FOLDS)
    inner_val, inner_train = np.arange(cut), np.arange(cut, len(targets))
    if len(inner_train) < 2:  # pragma: no cover - guarded by the size check
        return len(RELATIVE_PENALTIES) // 2
    average = targets[inner_train].mean(axis=0)
    candidates = _ridge_predictions(
        features[inner_train], targets[inner_train] - average, features[inner_val]
    )
    accuracy = [
        np.mean([_top_direction(average + c[i], targets[g]) for i, g in enumerate(inner_val)])
        for c in candidates
    ]
    return int(np.argmax(accuracy))


def benchmark_matrix(
    features: np.ndarray, targets: np.ndarray, *, folds: int = FOLDS, seed: int = 0
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-held-out-gene accuracy and correlation for the baseline and the model.

    Returns ``(baseline_accuracy, model_accuracy, baseline_corr, model_corr)``,
    one entry per gene.
    """
    if len(features) != len(targets):
        raise ValueError(f"{len(features)} features for {len(targets)} targets.")
    if len(targets) < folds * 2:
        raise ValueError(f"need at least {folds * 2} genes for {folds}-fold CV.")
    order = np.random.default_rng(seed).permutation(len(targets))
    base_acc, model_acc = np.empty(len(targets)), np.empty(len(targets))
    base_cor, model_cor = np.empty(len(targets)), np.empty(len(targets))

    for fold in range(folds):
        test = order[fold::folds]
        train = np.setdiff1d(order, test)
        average = targets[train].mean(axis=0)
        chosen = _choose_penalty(features[train], targets[train])
        # The model predicts the *departure* from the average, so any gain it
        # shows is gain over the baseline rather than a rediscovery of it.
        predicted = _ridge_predictions(features[train], targets[train] - average, features[test])[
            chosen
        ]
        for position, gene in enumerate(test):
            actual = targets[gene]
            model = average + predicted[position]
            base_acc[gene] = _top_direction(average, actual)
            base_cor[gene] = float(np.corrcoef(average, actual)[0, 1])
            model_acc[gene] = _top_direction(model, actual)
            model_cor[gene] = float(np.corrcoef(model, actual)[0, 1])
    return base_acc, model_acc, base_cor, model_cor


def trainable(screen: Screen) -> tuple[np.ndarray, np.ndarray]:
    """Features and targets for genes that were both knocked down and read out.

    A gene needs a row (what its knockdown did) to train on and a column (how it
    reacts to others) to be described by, so only genes with both can be used.
    """
    column = {gene: i for i, gene in enumerate(screen.readout_genes)}
    rows = [r for r in np.flatnonzero(screen.scorable) if screen.perturbed_genes[r] in column]
    if not rows:
        raise ValueError(f"{screen.key}: no gene is both knocked down and read out.")
    take = np.array([column[screen.perturbed_genes[r]] for r in rows])
    targets = screen.responses[np.array(rows)]
    features = screen.responses[:, take].T
    features = (features - features.mean(axis=0)) / (features.std(axis=0) + 1e-9)
    return features, targets


def run(screen: Screen, *, folds: int = FOLDS, seed: int = 0) -> Benchmark:
    features, targets = trainable(screen)
    base_acc, model_acc, base_cor, model_cor = benchmark_matrix(
        features, targets, folds=folds, seed=seed
    )
    difference = stats.ttest_rel(model_acc, base_acc)
    return Benchmark(
        screen=screen.key,
        cell_line=screen.cell_line,
        n_genes=len(targets),
        baseline_accuracy=round(float(base_acc.mean()), 4),
        model_accuracy=round(float(model_acc.mean()), 4),
        baseline_correlation=round(float(base_cor.mean()), 4),
        model_correlation=round(float(model_cor.mean()), 4),
        gain_pp=round(float((model_acc - base_acc).mean() * 100), 3),
        gain_p=float(difference.pvalue),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folds", type=int, default=FOLDS)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", type=Path, default=None, help="write the results as JSON")
    args = parser.parse_args(argv)
    results = []
    for key in SCREENS:
        screen = read_screen(fetch_artifact(key), key)
        results.append(asdict(run(screen, folds=args.folds, seed=args.seed)))
        print(json.dumps(results[-1], indent=2), flush=True)
    if args.out:
        args.out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
