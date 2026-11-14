"""Scores for held-out predictions of drug response.

r2                     share of the variation in response that is explained.
                       1 = perfect, 0 = no better than predicting the average,
                       negative = worse than that.
spearman               rank correlation: does the model order cell lines correctly?
within_group_spearman  the rank correlation inside each cancer type, then the
                       median over types. A model that only tells cancer types
                       apart scores near zero here.
bootstrap_gain         a 95% interval for the difference in R2 between two
                       models, from resampling cell lines.
repeated_gain          the same when each model has one set of held-out
                       predictions per repeat of the cross-validation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def r2(observed: np.ndarray, predicted: np.ndarray) -> float:
    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if observed.shape != predicted.shape:
        raise ValueError("observed and predicted must have the same shape")
    if np.isnan(observed).any() or np.isnan(predicted).any():
        raise ValueError("NaN in input")
    total = np.sum((observed - observed.mean()) ** 2)
    if total == 0:
        return float("nan")
    return float(1.0 - np.sum((observed - predicted) ** 2) / total)


def spearman(observed: np.ndarray, predicted: np.ndarray) -> float:
    if len(observed) < 3 or np.std(predicted) == 0 or np.std(observed) == 0:
        return float("nan")
    return float(spearmanr(observed, predicted).statistic)


def within_group_spearman(observed: pd.Series, predicted: pd.Series, groups: pd.Series, minimum: int = 20) -> float:
    """Median over groups (with at least `minimum` members) of the within-group Spearman."""
    values = []
    for _, index in groups.groupby(groups).groups.items():
        if len(index) >= minimum:
            values.append(spearman(observed.loc[index].to_numpy(), predicted.loc[index].to_numpy()))
    values = [v for v in values if not np.isnan(v)]
    return float(np.median(values)) if values else float("nan")


def bootstrap_gain(observed: np.ndarray, better: np.ndarray, simpler: np.ndarray,
                   n_resamples: int = 2000, seed: int = 0) -> tuple[float, float, float]:
    """R2(better) - R2(simpler), with a 95% interval from resampling cell lines."""
    observed, better, simpler = (np.asarray(a, dtype=float) for a in (observed, better, simpler))
    rng = np.random.default_rng(seed)
    n = len(observed)
    gains = np.empty(n_resamples)
    for i in range(n_resamples):
        pick = rng.integers(0, n, n)
        gains[i] = r2(observed[pick], better[pick]) - r2(observed[pick], simpler[pick])
    point = r2(observed, better) - r2(observed, simpler)
    low, high = np.nanpercentile(gains, [2.5, 97.5])
    return float(point), float(low), float(high)


def repeated_gain(observed: np.ndarray, better: np.ndarray, simpler: np.ndarray,
                  n_resamples: int = 2000, seed: int = 0) -> tuple[float, float, float]:
    """Mean over repeats of R2(better) - R2(simpler), with a 95% interval.

    `better` and `simpler` have one row per cell line and one column per
    repeat of the cross-validation. The gain is computed within each repeat
    and then averaged, so the point estimate equals the difference between the
    two mean R2 values in the results table. (Averaging the predictions over
    repeats first would flatter the less stable model, because an average of
    five predictions is a small ensemble.) The interval resamples cell lines,
    using the same resample for every repeat. It reflects which cell lines
    happened to be in the study, not how much a model changes if retrained.
    """
    observed = np.asarray(observed, dtype=float)
    better, simpler = np.asarray(better, dtype=float), np.asarray(simpler, dtype=float)
    if better.shape != simpler.shape or better.shape[0] != len(observed):
        raise ValueError("better and simpler must both be (cell lines x repeats)")

    def gain(rows: np.ndarray) -> float:
        return float(np.mean([r2(observed[rows], better[rows, j]) - r2(observed[rows], simpler[rows, j])
                              for j in range(better.shape[1])]))

    rng = np.random.default_rng(seed)
    n = len(observed)
    resampled = np.array([gain(rng.integers(0, n, n)) for _ in range(n_resamples)])
    low, high = np.nanpercentile(resampled, [2.5, 97.5])
    return gain(np.arange(n)), float(low), float(high)
