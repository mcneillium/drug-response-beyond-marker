"""Experiment 3: are groups found by clustering stable, and what drives them?

A clustering method always returns groups, whether or not any exist. Three
checks are made before believing them:

    stability   do the same groups come back when 20% of the cell lines are
                left out? (adjusted Rand index against the full clustering)
    reference   is that stability clearly higher than for data that has no
                groups at all?
    drivers     do the groups simply repeat tissue or sequencing centre?
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score

from drbm import models

N_GENES = 2000
N_COMPONENTS = 20
K_VALUES = (2, 3, 4, 5)
N_SUBSAMPLES = 100
SUBSAMPLE_SHARE = 0.80
N_REFERENCE = 20
STABLE, UNSTABLE = 0.75, 0.50


def effect_design(meta: pd.DataFrame, columns: list[str]) -> np.ndarray:
    """An intercept plus one 0/1 column per level of each named column (first level dropped)."""
    dummies = pd.get_dummies(meta[columns], drop_first=True, dtype=float)
    return np.column_stack([np.ones(len(meta)), dummies.to_numpy()])


def remove_effects(values: pd.DataFrame, meta: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Subtract the average effect of the named factors from every column of `values`.

    Fitted by least squares on all rows given. What is left (the residual) has
    the same mean in every tissue and in every sequencing centre.
    """
    design = effect_design(meta.loc[values.index], columns)
    matrix = values.to_numpy(dtype=float)
    coefficients, *_ = np.linalg.lstsq(design, matrix, rcond=None)
    return pd.DataFrame(matrix - design @ coefficients, index=values.index, columns=values.columns)


def components(expression: pd.DataFrame, n_genes: int = N_GENES, n_components: int = N_COMPONENTS) -> np.ndarray:
    """Most variable genes, scaled, reduced to principal components."""
    genes = models.top_variable(expression, min(n_genes, expression.shape[1]))
    scaled, _ = models.standardise(expression[genes], expression[genes])
    n_components = min(n_components, scaled.shape[0] - 1, scaled.shape[1])
    return PCA(n_components=n_components, random_state=0).fit_transform(scaled)


def kmeans(points: np.ndarray, k: int, seed: int = 0) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(points)


def stability(make_points, n_rows: int, k: int, n_subsamples: int = N_SUBSAMPLES,
              share: float = SUBSAMPLE_SHARE, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Cluster everything, then cluster random subsamples and compare.

    `make_points(row_positions)` returns the points to cluster for those rows,
    so that every preparation step is repeated on each subsample.

    Returns the full clustering and one adjusted Rand index per subsample.
    """
    everything = np.arange(n_rows)
    full = kmeans(make_points(everything), k, seed)
    rng = np.random.default_rng(seed)
    size = int(round(share * n_rows))
    agreement = np.empty(n_subsamples)
    for i in range(n_subsamples):
        rows = np.sort(rng.choice(n_rows, size=size, replace=False))
        agreement[i] = adjusted_rand_score(full[rows], kmeans(make_points(rows), k, seed + 1 + i))
    return full, agreement


def reference_stability(points: np.ndarray, k: int, n_reference: int = N_REFERENCE,
                        n_subsamples: int = 30, seed: int = 0) -> np.ndarray:
    """Stability of k-means on data with NO groups but the same spread per component.

    Each reference dataset is a single bell-shaped cloud: independent normal
    draws with the standard deviation of each column of `points`. Returns the
    mean adjusted Rand index of each reference dataset.
    """
    rng = np.random.default_rng(seed)
    spread = points.std(axis=0, ddof=1)
    means = np.empty(n_reference)
    for i in range(n_reference):
        cloud = rng.normal(0.0, spread, size=points.shape)
        _, agreement = stability(lambda rows, cloud=cloud: cloud[rows], len(cloud), k,
                                 n_subsamples=n_subsamples, seed=seed + 1000 * (i + 1))
        means[i] = agreement.mean()
    return means


def tissue_agreement(tissue: pd.Series, labels: np.ndarray, minimum: int = 5, n_shuffles: int = 200,
                     seed: int = 0) -> tuple[float, float]:
    """Share of cell lines that sit in the commonest group of their own tissue.

    Only tissues with at least `minimum` cell lines count (a tissue with one
    cell line would always score 100%). Returns the share, and the share
    expected if the groups had nothing to do with tissue (labels shuffled).

    The adjusted Rand index against tissue is low whenever a few groups each
    swallow several whole tissues. This measure is not: it is 1.0 in that case.
    """
    tissue = pd.Series(np.asarray(tissue))
    counts = tissue.value_counts()
    use = tissue.isin(counts[counts >= minimum].index).to_numpy()
    tissue, labels = tissue[use].to_numpy(), np.asarray(labels)[use]

    def share(groups: np.ndarray) -> float:
        return float(pd.crosstab(tissue, groups).max(axis=1).sum() / len(tissue))

    rng = np.random.default_rng(seed)
    shuffled = [share(rng.permutation(labels)) for _ in range(n_shuffles)]
    return share(labels), float(np.mean(shuffled))


def verdict(mean_agreement: float) -> str:
    if mean_agreement >= STABLE:
        return "stable"
    return "unstable" if mean_agreement < UNSTABLE else "in between"
