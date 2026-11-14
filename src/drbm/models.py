"""The baseline ladder: each step adds one kind of information.

    L0  nothing              the training-fold mean
    L1  PTEN status          two group means
    L2  background factors   ridge regression on cancer type and sequencing centre
    L3  PTEN + background    ridge regression
    L4  L3 + gene expression elastic net
    L5  L3 + gene expression gradient boosting

Every function here takes training rows and test rows and returns predictions
for the test rows. Anything learned from data (which genes vary most, how to
scale them, how strongly to regularise) is learned from the training rows only.
"""

from __future__ import annotations

import warnings

import lightgbm as lgb
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV, RidgeCV

N_GENES = 2000
RIDGE_ALPHAS = np.logspace(-2, 4, 13)
GBM_PARAMS = {
    "objective": "regression",
    "n_estimators": 400,
    "learning_rate": 0.03,
    "num_leaves": 7,
    "min_child_samples": 10,
    "subsample": 0.8,
    "subsample_freq": 1,
    "colsample_bytree": 0.3,
    "reg_lambda": 1.0,
    "n_jobs": 2,
    "verbose": -1,
}
PATHWAY_GENES = ["PTEN", "AKT1", "AKT2", "AKT3", "PIK3CA", "PIK3CB", "PIK3R1", "MTOR", "TSC1", "TSC2", "INPP4B", "PDPK1"]
"""Declared in the protocol as a sanity check on the explanations."""


# ---------- building blocks ----------

def top_variable(frame: pd.DataFrame, n: int) -> list[str]:
    """The `n` columns with the largest variance. Does not look at the response."""
    variance = frame.var(axis=0, ddof=0)
    return variance.nlargest(n).index.tolist()


def standardise(train: pd.DataFrame, test: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Scale columns to mean 0 and spread 1 using the training rows only."""
    mean = train.mean(axis=0)
    spread = train.std(axis=0, ddof=0).replace(0, 1.0)
    return ((train - mean) / spread).to_numpy(dtype=float), ((test - mean) / spread).to_numpy(dtype=float)


def background(meta_train: pd.DataFrame, meta_test: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One column per cancer group and per sequencing centre seen in training."""
    def encode(meta: pd.DataFrame, columns: list[str] | None) -> pd.DataFrame:
        table = pd.get_dummies(meta[["cancer_group", "source"]], prefix=["type", "centre"], dtype=float)
        return table if columns is None else table.reindex(columns=columns, fill_value=0.0)
    train = encode(meta_train, None)
    return train, encode(meta_test, train.columns.tolist())


def pten_column(meta: pd.DataFrame, definition: str = "pten_altered") -> pd.DataFrame:
    return meta[[definition]].astype(float).rename(columns={definition: "PTEN_status"})


def ridge(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray) -> tuple[np.ndarray, float]:
    model = RidgeCV(alphas=RIDGE_ALPHAS).fit(train_x, train_y)
    return model.predict(test_x), float(model.alpha_)


# ---------- the ladder ----------

def l0_mean(y_train: pd.Series, test_ids: pd.Index) -> pd.Series:
    return pd.Series(y_train.mean(), index=test_ids)


def l1_pten(y_train: pd.Series, meta_train: pd.DataFrame, meta_test: pd.DataFrame,
            definition: str = "pten_altered") -> pd.Series:
    """The mean response of training cell lines with the same PTEN status."""
    means = y_train.groupby(meta_train[definition]).mean()
    return meta_test[definition].map(means).fillna(y_train.mean())


def l2_background(y_train, meta_train, meta_test) -> pd.Series:
    train, test = background(meta_train, meta_test)
    train_x, test_x = standardise(train, test)
    prediction, _ = ridge(train_x, y_train.to_numpy(), test_x)
    return pd.Series(prediction, index=meta_test.index)


def l3_pten_background(y_train, meta_train, meta_test, definition: str = "pten_altered") -> pd.Series:
    train_bg, test_bg = background(meta_train, meta_test)
    train = pd.concat([pten_column(meta_train, definition), train_bg], axis=1)
    test = pd.concat([pten_column(meta_test, definition), test_bg], axis=1)
    train_x, test_x = standardise(train, test)
    prediction, _ = ridge(train_x, y_train.to_numpy(), test_x)
    return pd.Series(prediction, index=meta_test.index)


def design_with_expression(expr_train, expr_test, meta_train, meta_test, n_genes: int = N_GENES):
    """PTEN + background + the most variable genes, as aligned training and test tables."""
    genes = top_variable(expr_train, n_genes)
    train_bg, test_bg = background(meta_train, meta_test)
    train = pd.concat([pten_column(meta_train), train_bg, expr_train[genes]], axis=1)
    test = pd.concat([pten_column(meta_test), test_bg, expr_test[genes]], axis=1)
    return train, test


def l4_elastic_net(y_train, expr_train, expr_test, meta_train, meta_test, seed: int = 0):
    """Elastic net on PTEN + background + genes. Returns predictions and non-zero weights."""
    train, test = design_with_expression(expr_train, expr_test, meta_train, meta_test)
    train_x, test_x = standardise(train, test)
    # 30 strengths from the largest useful value down to a hundredth of it. Random
    # coordinate order and a looser tolerance make each fit about 25 times faster
    # and gave the same held-out score in a timing check on one fold.
    model = ElasticNetCV(l1_ratio=0.5, alphas=30, eps=1e-2, cv=5, max_iter=5000, tol=1e-3,
                         selection="random", n_jobs=2, random_state=seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model.fit(train_x, y_train.to_numpy())
    weights = pd.Series(model.coef_, index=train.columns)
    return pd.Series(model.predict(test_x), index=meta_test.index), weights[weights != 0]


def l5_boosting(y_train, expr_train, expr_test, meta_train, meta_test, seed: int = 0):
    """Gradient boosting on the same inputs. Returns predictions and mean |SHAP| per feature."""
    train, test = design_with_expression(expr_train, expr_test, meta_train, meta_test)
    model = lgb.LGBMRegressor(**GBM_PARAMS, random_state=seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model.fit(train.to_numpy(dtype=float), y_train.to_numpy())
        prediction = model.predict(test.to_numpy(dtype=float))
        contributions = model.predict(test.to_numpy(dtype=float), pred_contrib=True)[:, :-1]
    importance = pd.Series(np.abs(contributions).mean(axis=0), index=train.columns)
    return pd.Series(prediction, index=meta_test.index), importance


# ---------- the leak demonstration ----------

def top_correlated(frame: pd.DataFrame, y: pd.Series, n: int) -> list[str]:
    """The `n` columns most correlated (in absolute value) with the response."""
    centred = frame - frame.mean(axis=0)
    target = (y - y.mean()).to_numpy()
    spread = np.sqrt((centred ** 2).sum(axis=0)) * np.sqrt((target ** 2).sum())
    correlation = (centred.to_numpy().T @ target) / spread.replace(0, np.nan).to_numpy()
    return pd.Series(np.abs(correlation), index=frame.columns).nlargest(n).index.tolist()


def ridge_on_genes(y_train, expr_train, expr_test, genes: list[str]) -> pd.Series:
    train_x, test_x = standardise(expr_train[genes], expr_test[genes])
    prediction, _ = ridge(train_x, y_train.to_numpy(), test_x)
    return pd.Series(prediction, index=expr_test.index)
