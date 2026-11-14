"""Experiment 7: how well do two measurements of the same drug response agree?

If a response is measured twice and each measurement is "truth + its own
noise", with equal and independent noise, then the Pearson correlation between
the two equals the share of the variation that is truth. No model predicting
one measurement can explain more than that share, however good it is. This is
called a noise ceiling.

The estimate is too low if the two measurements differ in more than noise (a
different dose range or a different assay), and too high if both share an
error (the same batch of cells, say).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


def pair_agreement(first: pd.Series, second: pd.Series, n_resamples: int = 0, seed: int = 0) -> dict:
    """Agreement between two measurements, on the cell lines that have both.

    Both series are indexed by cell line. With `n_resamples` > 0 a 95%
    interval for the Pearson correlation is added, from resampling cell lines.
    """
    shared = first.dropna().index.intersection(second.dropna().index)
    a, b = first.loc[shared].to_numpy(dtype=float), second.loc[shared].to_numpy(dtype=float)
    out = {"n_shared": len(shared), "pearson": float("nan"), "spearman": float("nan")}
    if len(shared) < 20 or a.std() == 0 or b.std() == 0:
        return out
    out["pearson"] = float(pearsonr(a, b).statistic)
    out["spearman"] = float(spearmanr(a, b).statistic)
    if n_resamples:
        rng = np.random.default_rng(seed)
        draws = []
        for _ in range(n_resamples):
            pick = rng.integers(0, len(a), len(a))
            draws.append(np.corrcoef(a[pick], b[pick])[0, 1])
        out["pearson_low"], out["pearson_high"] = (float(v) for v in np.nanpercentile(draws, [2.5, 97.5]))
    return out


def share_of_ceiling(model_r2: float, ceiling: float) -> float:
    """What fraction of the explainable variation a model's R2 amounts to."""
    if not ceiling > 0:
        return float("nan")
    return float(model_r2 / ceiling)
