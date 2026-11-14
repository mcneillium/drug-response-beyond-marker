import numpy as np
import pandas as pd
import pytest

from drbm import overlap, scores


def pieces(toy, hidden_test=(), hidden_train=()):
    ids = toy["ids"]
    train, test = ids[30:], ids[:30]
    avail = pd.Series(True, index=ids)
    avail[list(hidden_test) + list(hidden_train)] = False
    protein = toy["protein"].copy()
    return dict(y_train=toy["y"][train], expr_train=toy["expression"].loc[train], expr_test=toy["expression"].loc[test],
                prot_train=protein.loc[train], prot_test=protein.loc[test],
                avail_train=avail[train], avail_test=avail[test]), train, test


def test_keep_protein_hits_the_share_and_never_invents_data():
    has = pd.Series([True] * 90 + [False] * 10, index=[f"M{i}" for i in range(100)])
    kept = overlap.keep_protein(has, 0.5, seed=0)
    assert kept.sum() == 50
    assert not kept[~has].any()                         # lines without protein stay without it
    assert overlap.keep_protein(has, 0.95, seed=0).equals(has)   # cannot keep more than exists
    assert kept.equals(overlap.keep_protein(has, 0.5, seed=0))
    assert not kept.equals(overlap.keep_protein(has, 0.5, seed=1))
    with pytest.raises(ValueError):
        overlap.keep_protein(has, 0.0, seed=0)


def test_protein_tables_filter_and_fill_from_training_rows():
    train = pd.DataFrame({"often": [1.0, 2.0, 3.0, 4.0, np.nan], "rare": [1.0, np.nan, np.nan, np.nan, np.nan]})
    test = pd.DataFrame({"often": [np.nan], "rare": [9.0]})
    train_x, test_x, proteins = overlap.protein_tables(train, test)
    assert proteins == ["often"]                        # 'rare' is seen in 20% of training rows
    assert np.isfinite(train_x).all() and np.isfinite(test_x).all()
    # the missing test value was filled with the training median (2.5), which is also the mean
    assert test_x[0, 0] == pytest.approx(0.0)


def test_complete_cases_predicts_only_rows_with_protein(toy):
    ids = toy["ids"]
    args, train, test = pieces(toy, hidden_test=ids[:10], hidden_train=ids[30:60])
    prediction, _, n_used = overlap.complete_cases(**args)
    assert prediction[ids[:10]].isna().all()
    assert prediction[ids[10:30]].notna().all()
    assert n_used == 60


def test_complete_cases_never_reads_hidden_protein(toy):
    """Corrupting the protein values of hidden rows must change nothing."""
    ids = toy["ids"]
    args, _, _ = pieces(toy, hidden_test=ids[:10], hidden_train=ids[30:60])
    first, _, _ = overlap.complete_cases(**args)
    args["prot_train"] = args["prot_train"].copy()
    args["prot_train"].loc[ids[30:60]] = 1e6
    second, _, _ = overlap.complete_cases(**args)
    pd.testing.assert_series_equal(first, second)


def test_fill_gaps_predicts_everyone_and_never_reads_hidden_protein(toy):
    ids = toy["ids"]
    args, _, test = pieces(toy, hidden_test=ids[:10], hidden_train=ids[30:60])
    first, _ = overlap.fill_gaps(**args)
    assert first.notna().all() and first.index.equals(test)
    args["prot_train"], args["prot_test"] = args["prot_train"].copy(), args["prot_test"].copy()
    args["prot_train"].loc[ids[30:60]] = 1e6
    args["prot_test"].loc[ids[:10]] = -1e6
    second, _ = overlap.fill_gaps(**args)
    pd.testing.assert_series_equal(first, second)


def test_late_fusion_falls_back_to_expression_without_protein(toy):
    ids = toy["ids"]
    args, _, _ = pieces(toy, hidden_test=ids[:10], hidden_train=ids[30:60])
    fused, weight = overlap.late_fusion(**args)
    alone, _ = overlap.expression_only(args["y_train"], args["expr_train"], args["expr_test"])
    assert fused.notna().all()
    np.testing.assert_allclose(fused[ids[:10]], alone[ids[:10]])
    assert 0.0 <= weight <= 1.0


def test_late_fusion_weights_protein_when_protein_holds_the_signal(toy):
    """If the response is driven by a protein and not by any gene, the blend should lean on protein."""
    ids = toy["ids"]
    rng = np.random.default_rng(5)
    y = pd.Series(0.5 + 0.1 * toy["protein"]["p0"] + 0.01 * rng.normal(size=len(ids)), index=ids)
    args, train, test = pieces(toy)
    args["y_train"] = y[train]
    fused, weight = overlap.late_fusion(**args)
    alone, _ = overlap.expression_only(args["y_train"], args["expr_train"], args["expr_test"])
    assert weight >= 0.8
    assert scores.r2(y[test].to_numpy(), fused.to_numpy()) > scores.r2(y[test].to_numpy(), alone.to_numpy()) + 0.3


def test_choose_weight_extremes():
    y = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    noise = pd.Series(np.zeros(12))
    assert overlap.choose_weight(y, noise, y) == 1.0     # protein model is perfect
    assert overlap.choose_weight(y, y, noise) == 0.0     # expression model is perfect
    assert overlap.choose_weight(y, y, pd.Series(np.nan, index=y.index)) == 0.0   # no protein predictions at all
