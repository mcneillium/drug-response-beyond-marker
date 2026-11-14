"""Experiment 2: what to do when one data type is missing for some cell lines.

Gene expression is available for every cell line. Protein data is available
for some. Four strategies, all using ridge regression so only the strategy
differs:

    expression_only   ignore protein altogether
    complete_cases    train and predict only on cell lines that have both
    fill_gaps         give cell lines without protein the training median of
                      each protein, plus a flag saying "this row was filled"
    late_fusion       one model per data type, each trained on the cell lines
                      that have it; blend the two predictions where both exist

Everything learned from data (genes, proteins, medians, scaling, ridge
strength, blend weight) is learned from the training rows only.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

from drbm import models

ALPHAS = np.logspace(0, 6, 13)
MIN_DETECTED = 0.80                     # a protein must be seen in this share of training lines with protein
BLEND_WEIGHTS = np.round(np.arange(0.0, 1.01, 0.1), 2)
STRATEGIES = ("expression_only", "complete_cases", "fill_gaps", "late_fusion")


def keep_protein(has_protein: pd.Series, share: float, seed: int) -> pd.Series:
    """Choose which cell lines keep their protein data.

    `share` is the share of ALL cell lines that should have protein afterwards.
    Cell lines that never had protein stay without it. If the share asked for
    is at or above what exists, nothing is hidden.
    """
    if not 0 < share <= 1:
        raise ValueError("share must be in (0, 1]")
    have = has_protein[has_protein].index
    target = int(round(share * len(has_protein)))
    if target >= len(have):
        return has_protein.copy()
    rng = np.random.default_rng(seed)
    kept = rng.choice(have.to_numpy(), size=target, replace=False)
    return pd.Series(has_protein.index.isin(kept), index=has_protein.index)


def _ridge(train_x: np.ndarray, y: np.ndarray, test_x: np.ndarray) -> tuple[np.ndarray, float]:
    model = RidgeCV(alphas=ALPHAS).fit(train_x, y)
    return model.predict(test_x), float(model.alpha_)


def _scale(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean 0 and spread 1 per column, from the training rows only. Plain arrays, for speed."""
    mean = train.mean(axis=0)
    spread = train.std(axis=0)
    spread[spread == 0] = 1.0
    return (train - mean) / spread, (test - mean) / spread


def _expression_tables(expr_train: pd.DataFrame, expr_test: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    genes = models.top_variable(expr_train, models.N_GENES)
    return _scale(expr_train[genes].to_numpy(dtype=float), expr_test[genes].to_numpy(dtype=float))


def _protein_arrays(prot_train: pd.DataFrame, prot_test: pd.DataFrame, reference_rows: np.ndarray | None = None):
    """Filter proteins and fill gaps, learning both from the training rows that have protein data.

    `reference_rows` marks those rows (default: all training rows). Returns the
    filled training and test arrays, unscaled, and the names of the proteins kept.
    """
    train = prot_train.to_numpy(dtype=float)
    test = prot_test.to_numpy(dtype=float)
    reference = train if reference_rows is None else train[reference_rows]
    keep = (~np.isnan(reference)).mean(axis=0) >= MIN_DETECTED
    train, test, reference = train[:, keep], test[:, keep], reference[:, keep]
    median = np.nanmedian(reference, axis=0)
    train = np.where(np.isnan(train), median, train)
    test = np.where(np.isnan(test), median, test)
    return train, test, prot_train.columns[keep].tolist()


def protein_tables(prot_train: pd.DataFrame, prot_test: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Protein inputs for rows that have protein data.

    Keeps proteins detected in at least 80% of the training rows, fills the
    remaining gaps with the training median, and scales with training values.
    """
    train, test, proteins = _protein_arrays(prot_train, prot_test)
    train_x, test_x = _scale(train, test)
    return train_x, test_x, proteins


def expression_only(y_train, expr_train, expr_test) -> tuple[pd.Series, float]:
    train_x, test_x = _expression_tables(expr_train, expr_test)
    prediction, alpha = _ridge(train_x, y_train.to_numpy(), test_x)
    return pd.Series(prediction, index=expr_test.index), alpha


def complete_cases(y_train, expr_train, expr_test, prot_train, prot_test, avail_train, avail_test):
    """Train on training rows with both data types; predict test rows with both. Others get NaN."""
    rows, targets = avail_train[avail_train].index, avail_test[avail_test].index
    out = pd.Series(np.nan, index=expr_test.index)
    if len(targets) == 0:
        return out, float("nan"), len(rows)
    expr_x, expr_t = _expression_tables(expr_train.loc[rows], expr_test.loc[targets])
    prot_x, prot_t, _ = protein_tables(prot_train.loc[rows], prot_test.loc[targets])
    prediction, alpha = _ridge(np.hstack([expr_x, prot_x]), y_train[rows].to_numpy(), np.hstack([expr_t, prot_t]))
    out.loc[targets] = prediction
    return out, alpha, len(rows)


def fill_gaps(y_train, expr_train, expr_test, prot_train, prot_test, avail_train, avail_test):
    """Train on every training row; rows without protein get training medians and a flag."""
    have_train, have_test = avail_train.to_numpy(dtype=bool), avail_test.to_numpy(dtype=bool)
    # Hide the protein values of rows that are meant to have none, so they cannot be read.
    hidden_train = prot_train.where(avail_train, np.nan)
    hidden_test = prot_test.where(avail_test, np.nan)
    train, test, _ = _protein_arrays(hidden_train, hidden_test, reference_rows=have_train)
    train = np.column_stack([train, (~have_train).astype(float)])
    test = np.column_stack([test, (~have_test).astype(float)])
    prot_x, prot_t = _scale(train, test)
    expr_x, expr_t = _expression_tables(expr_train, expr_test)
    prediction, alpha = _ridge(np.hstack([expr_x, prot_x]), y_train.to_numpy(), np.hstack([expr_t, prot_t]))
    return pd.Series(prediction, index=expr_test.index), alpha


def _protein_model(y_train, prot_train, prot_test, avail_train, avail_test) -> pd.Series:
    """Predictions from a model trained on the rows with protein. NaN for rows without."""
    out = pd.Series(np.nan, index=prot_test.index)
    rows, targets = avail_train[avail_train].index, avail_test[avail_test].index
    if len(rows) >= 20 and len(targets) > 0:
        prot_x, prot_t, _ = protein_tables(prot_train.loc[rows], prot_test.loc[targets])
        out.loc[targets] = _ridge(prot_x, y_train[rows].to_numpy(), prot_t)[0]
    return out


def _inner_folds(y_train: pd.Series, seed: int):
    splitter = KFold(n_splits=5, shuffle=True, random_state=seed)
    for fit_positions, hold_positions in splitter.split(np.zeros(len(y_train))):
        yield y_train.index[fit_positions], y_train.index[hold_positions]


def expression_part(y_train, expr_train, expr_test, seed: int = 0) -> tuple[pd.Series, pd.Series]:
    """The expression half of late fusion: inner held-out predictions, and test predictions.

    It does not depend on which cell lines have protein, so it can be computed
    once per fold and reused for every share.
    """
    inner = pd.Series(np.nan, index=y_train.index)
    for fit, hold in _inner_folds(y_train, seed):
        inner.loc[hold] = expression_only(y_train[fit], expr_train.loc[fit], expr_train.loc[hold])[0]
    return inner, expression_only(y_train, expr_train, expr_test)[0]


def choose_weight(y: pd.Series, from_expression: pd.Series, from_protein: pd.Series) -> float:
    """The blend weight on the protein model with the smallest squared error."""
    both = from_protein.notna()
    if both.sum() < 10:
        return 0.0
    errors = [float(np.mean((y[both] - (w * from_protein[both] + (1 - w) * from_expression[both])) ** 2))
              for w in BLEND_WEIGHTS]
    return float(BLEND_WEIGHTS[int(np.argmin(errors))])


def late_fusion(y_train, expr_train, expr_test, prot_train, prot_test, avail_train, avail_test, seed: int = 0,
                expression=None):
    """Blend an expression model and a protein model. Returns predictions and the weight on protein.

    The weight is chosen on held-out predictions made inside the training
    rows. `expression` may carry a precomputed `expression_part` for the same
    rows and seed.
    """
    inner_expression, from_expression = expression or expression_part(y_train, expr_train, expr_test, seed)
    inner_protein = pd.Series(np.nan, index=y_train.index)
    for fit, hold in _inner_folds(y_train, seed):
        inner_protein.loc[hold] = _protein_model(y_train[fit], prot_train.loc[fit], prot_train.loc[hold],
                                                 avail_train[fit], avail_train[hold])
    weight = choose_weight(y_train, inner_expression, inner_protein)
    from_protein = _protein_model(y_train, prot_train, prot_test, avail_train, avail_test)
    blended = weight * from_protein + (1 - weight) * from_expression
    return blended.where(from_protein.notna(), from_expression), weight
