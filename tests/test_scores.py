import numpy as np
import pandas as pd
import pytest

from drbm import scores


def test_r2_hand_worked():
    observed = np.array([1.0, 2.0, 3.0, 4.0])          # mean 2.5, total sum of squares 5
    predicted = np.array([1.0, 2.0, 3.0, 2.0])         # one error of 2 -> squared error 4
    assert scores.r2(observed, predicted) == pytest.approx(1 - 4 / 5)


def test_r2_reference_points():
    observed = np.array([0.2, 0.5, 0.9, 0.4])
    assert scores.r2(observed, observed) == pytest.approx(1.0)
    assert scores.r2(observed, np.full(4, observed.mean())) == pytest.approx(0.0)
    assert scores.r2(observed, np.full(4, 5.0)) < 0


def test_r2_refuses_missing_values_and_shape_mismatch():
    with pytest.raises(ValueError):
        scores.r2(np.array([1.0, np.nan]), np.array([1.0, 2.0]))
    with pytest.raises(ValueError):
        scores.r2(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0]))


def test_r2_of_constant_observations_is_undefined():
    assert np.isnan(scores.r2(np.ones(5), np.arange(5.0)))


def test_spearman_only_cares_about_order():
    observed = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert scores.spearman(observed, observed ** 3) == pytest.approx(1.0)
    assert scores.spearman(observed, -observed) == pytest.approx(-1.0)
    assert np.isnan(scores.spearman(observed, np.ones(5)))


def test_within_group_spearman_ignores_differences_between_groups():
    """A model that only knows the group gets the overall order partly right but nothing within groups."""
    rng = np.random.default_rng(1)
    groups = pd.Series(np.repeat(["a", "b"], 50))
    observed = pd.Series(np.where(groups == "a", 0.0, 10.0) + rng.normal(size=100))
    knows_group_only = pd.Series(np.where(groups == "a", 0.0, 10.0) + 1e-6 * rng.normal(size=100))
    assert scores.spearman(observed.to_numpy(), knows_group_only.to_numpy()) > 0.7
    assert abs(scores.within_group_spearman(observed, knows_group_only, groups)) < 0.3
    assert scores.within_group_spearman(observed, observed, groups) == pytest.approx(1.0)


def test_within_group_spearman_skips_small_groups():
    groups = pd.Series(["big"] * 30 + ["small"] * 5)
    observed = pd.Series(np.arange(35.0))
    predicted = pd.Series(np.r_[np.arange(30.0), -np.arange(5.0)])   # small group is predicted backwards
    assert scores.within_group_spearman(observed, predicted, groups, minimum=20) == pytest.approx(1.0)


def test_bootstrap_gain_point_and_interval():
    rng = np.random.default_rng(2)
    observed = rng.normal(size=300)
    better = observed + 0.3 * rng.normal(size=300)
    simpler = np.zeros(300) + observed.mean()
    point, low, high = scores.bootstrap_gain(observed, better, simpler, n_resamples=300)
    assert point == pytest.approx(scores.r2(observed, better) - scores.r2(observed, simpler))
    assert low < point < high and low > 0
    again = scores.bootstrap_gain(observed, better, simpler, n_resamples=300)
    assert again == (point, low, high)                   # same seed, same answer


def test_bootstrap_gain_of_identical_models_is_zero():
    rng = np.random.default_rng(3)
    observed = rng.normal(size=100)
    prediction = observed + rng.normal(size=100)
    assert scores.bootstrap_gain(observed, prediction, prediction, n_resamples=100) == (0.0, 0.0, 0.0)


def test_repeated_gain_matches_the_difference_of_mean_r2():
    rng = np.random.default_rng(4)
    observed = rng.normal(size=200)
    better = observed[:, None] + 0.5 * rng.normal(size=(200, 5))
    simpler = observed[:, None] + 1.0 * rng.normal(size=(200, 5))
    point, low, high = scores.repeated_gain(observed, better, simpler, n_resamples=200)
    expected = np.mean([scores.r2(observed, better[:, j]) for j in range(5)]) - \
        np.mean([scores.r2(observed, simpler[:, j]) for j in range(5)])
    assert point == pytest.approx(expected)
    assert low < point < high and low > 0


def test_repeated_gain_does_not_reward_averaging_noisy_predictions():
    """Two equally noisy models: averaging one's repeats first would make it look better. This must not."""
    rng = np.random.default_rng(5)
    observed = rng.normal(size=400)
    noisy = observed[:, None] + rng.normal(size=(400, 5))            # differs from repeat to repeat
    steady = np.repeat((observed + rng.normal(size=400))[:, None], 5, axis=1)   # same in every repeat
    point, low, high = scores.repeated_gain(observed, noisy, steady, n_resamples=200)
    assert low < 0 < high
    flattered = scores.r2(observed, noisy.mean(axis=1)) - scores.r2(observed, steady.mean(axis=1))
    assert flattered > 0.5                                # the method this function replaces


def test_repeated_gain_checks_shapes():
    with pytest.raises(ValueError):
        scores.repeated_gain(np.zeros(5), np.zeros((5, 2)), np.zeros((5, 3)))
