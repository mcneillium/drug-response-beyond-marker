import numpy as np
import pandas as pd
import pytest

from drbm import repeats


def two_measurements(n=4000, noise=0.5, seed=0):
    rng = np.random.default_rng(seed)
    ids = [f"M{i}" for i in range(n)]
    truth = rng.normal(size=n)
    return (pd.Series(truth + noise * rng.normal(size=n), index=ids),
            pd.Series(truth + noise * rng.normal(size=n), index=ids), truth)


def test_correlation_of_repeats_equals_the_share_that_is_real():
    first, second, _ = two_measurements(noise=0.5)
    expected = 1.0 / (1.0 + 0.5 ** 2)                    # real variance / total variance = 0.8
    assert repeats.pair_agreement(first, second)["pearson"] == pytest.approx(expected, abs=0.02)


def test_the_ceiling_really_is_a_ceiling():
    """Even a model that knows the truth exactly cannot beat the ceiling against a noisy measurement."""
    from drbm import scores
    first, second, truth = two_measurements(noise=0.5)
    ceiling = repeats.pair_agreement(first, second)["pearson"]
    assert scores.r2(first.to_numpy(), truth) == pytest.approx(ceiling, abs=0.03)


def test_only_shared_cell_lines_are_used():
    first = pd.Series(np.arange(30.0), index=[f"M{i}" for i in range(30)])
    second = pd.Series(np.arange(40.0), index=[f"M{i}" for i in range(10, 50)])
    second.iloc[0] = np.nan
    result = repeats.pair_agreement(first, second)
    assert result["n_shared"] == 19                      # M11..M29; M10 is missing in the second
    assert np.isnan(result["pearson"])                   # fewer than 20: no estimate


def test_interval_brackets_the_estimate_and_is_reproducible():
    first, second, _ = two_measurements(n=300)
    result = repeats.pair_agreement(first, second, n_resamples=200)
    assert result["pearson_low"] < result["pearson"] < result["pearson_high"]
    assert result == repeats.pair_agreement(first, second, n_resamples=200)


def test_share_of_ceiling():
    assert repeats.share_of_ceiling(0.2, 0.8) == pytest.approx(0.25)
    assert np.isnan(repeats.share_of_ceiling(0.2, 0.0))
