import numpy as np
import pandas as pd
import pytest

from drbm import models, scores


def split(toy, n_test=30):
    ids = toy["ids"]
    return ids[n_test:], ids[:n_test]


def test_top_variable_picks_by_variance_and_ignores_the_response():
    frame = pd.DataFrame({"flat": [1.0, 1.0, 1.0, 1.0], "wide": [0.0, 10.0, 0.0, 10.0], "mid": [0.0, 1.0, 0.0, 1.0]})
    assert models.top_variable(frame, 2) == ["wide", "mid"]


def test_standardise_uses_training_rows_only():
    train = pd.DataFrame({"a": [0.0, 2.0], "constant": [5.0, 5.0]})
    test = pd.DataFrame({"a": [100.0], "constant": [5.0]})
    train_x, test_x = models.standardise(train, test)
    assert train_x[:, 0].tolist() == [-1.0, 1.0]
    assert test_x[0, 0] == pytest.approx(99.0)          # scaled with the training mean 1 and spread 1
    assert np.isfinite(train_x).all() and np.isfinite(test_x).all()   # constant column does not divide by zero


def test_background_columns_come_from_training_rows(toy):
    train, test = split(toy)
    meta_test = toy["meta"].loc[test].copy()
    meta_test.loc[test[0], "cancer_group"] = "Never seen in training"
    train_table, test_table = models.background(toy["meta"].loc[train], meta_test)
    assert train_table.columns.tolist() == test_table.columns.tolist()
    assert test_table.loc[test[0]].filter(like="type_").sum() == 0     # unseen type gets no column


def test_l0_and_l1_are_training_means(toy):
    train, test = split(toy)
    y, meta = toy["y"], toy["meta"]
    assert (models.l0_mean(y[train], test) == y[train].mean()).all()
    prediction = models.l1_pten(y[train], meta.loc[train], meta.loc[test])
    altered_mean = y[train][meta.loc[train, "pten_altered"]].mean()
    normal_mean = y[train][~meta.loc[train, "pten_altered"]].mean()
    assert prediction[meta.loc[test, "pten_altered"]].eq(altered_mean).all()
    assert prediction[~meta.loc[test, "pten_altered"]].eq(normal_mean).all()


def test_l1_recovers_a_real_marker_effect(toy):
    train, test = split(toy)
    meta = toy["meta"]
    y = pd.Series(np.where(meta["pten_altered"], 0.4, 0.9), index=toy["ids"])
    prediction = models.l1_pten(y[train], meta.loc[train], meta.loc[test])
    assert scores.r2(y[test].to_numpy(), prediction.to_numpy()) == pytest.approx(1.0)


def test_ladder_steps_return_one_prediction_per_test_row(toy):
    train, test = split(toy)
    y, meta, expr = toy["y"], toy["meta"], toy["expression"]
    for prediction in (
        models.l2_background(y[train], meta.loc[train], meta.loc[test]),
        models.l3_pten_background(y[train], meta.loc[train], meta.loc[test]),
        models.l4_elastic_net(y[train], expr.loc[train], expr.loc[test], meta.loc[train], meta.loc[test])[0],
        models.l5_boosting(y[train], expr.loc[train], expr.loc[test], meta.loc[train], meta.loc[test])[0],
    ):
        assert prediction.index.equals(test)
        assert prediction.notna().all()


def test_expression_models_find_the_gene_that_drives_the_response(toy):
    train, test = split(toy)
    y, meta, expr = toy["y"], toy["meta"], toy["expression"]
    prediction, weights = models.l4_elastic_net(y[train], expr.loc[train], expr.loc[test], meta.loc[train], meta.loc[test])
    assert scores.r2(y[test].to_numpy(), prediction.to_numpy()) > 0.8
    assert weights.abs().idxmax() == "g0"
    prediction, importance = models.l5_boosting(y[train], expr.loc[train], expr.loc[test], meta.loc[train], meta.loc[test])
    assert scores.r2(y[test].to_numpy(), prediction.to_numpy()) > 0.5
    assert importance.idxmax() == "g0"


def test_design_with_expression_chooses_genes_from_training_rows_only(toy):
    """Changing the test rows must not change which genes are used."""
    train, test = split(toy)
    meta, expr = toy["meta"], toy["expression"]
    first, _ = models.design_with_expression(expr.loc[train], expr.loc[test], meta.loc[train], meta.loc[test], n_genes=10)
    altered_test = expr.loc[test] * 1000.0
    second, _ = models.design_with_expression(expr.loc[train], altered_test, meta.loc[train], meta.loc[test], n_genes=10)
    assert first.columns.tolist() == second.columns.tolist()


def test_top_correlated_matches_numpy(toy):
    expr, y = toy["expression"], toy["y"]
    chosen = models.top_correlated(expr, y, 5)
    correlations = expr.apply(lambda column: abs(np.corrcoef(column, y)[0, 1]))
    assert chosen == correlations.nlargest(5).index.tolist()
    assert chosen[0] == "g0"


def test_the_leak_on_pure_noise():
    """With nothing to find, choosing genes on all rows first still looks like skill. Choosing in-fold does not."""
    rng = np.random.default_rng(0)
    ids = pd.Index([f"M{i}" for i in range(100)])
    expr = pd.DataFrame(rng.normal(size=(100, 3000)), index=ids, columns=[f"g{j}" for j in range(3000)])
    y = pd.Series(rng.normal(size=100), index=ids)
    leaked_genes = models.top_correlated(expr, y, 50)
    wrong, right = [], []
    for start in range(0, 100, 20):
        test = ids[start:start + 20]
        train = ids.difference(test)
        wrong.append(models.ridge_on_genes(y[train], expr.loc[train], expr.loc[test], leaked_genes))
        honest_genes = models.top_correlated(expr.loc[train], y[train], 50)
        right.append(models.ridge_on_genes(y[train], expr.loc[train], expr.loc[test], honest_genes))
    wrong, right = pd.concat(wrong), pd.concat(right)
    assert scores.r2(y[wrong.index].to_numpy(), wrong.to_numpy()) > 0.2
    assert scores.r2(y[right.index].to_numpy(), right.to_numpy()) < 0.05
