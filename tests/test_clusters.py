import numpy as np
import pandas as pd
import pytest

from drbm import clusters


def blobs(n_per_group=40, n_genes=50, separation=8.0, seed=0):
    rng = np.random.default_rng(seed)
    centres = np.zeros((2, n_genes))
    centres[1, :10] = separation
    values = np.vstack([centres[g] + rng.normal(size=(n_per_group, n_genes)) for g in (0, 1)])
    index = pd.Index([f"M{i}" for i in range(2 * n_per_group)])
    return pd.DataFrame(values, index=index, columns=[f"g{j}" for j in range(n_genes)]), np.repeat([0, 1], n_per_group)


def test_remove_effects_makes_group_means_equal():
    rng = np.random.default_rng(0)
    meta = pd.DataFrame({"tissue": np.repeat(["a", "b", "c"], 30), "source": np.tile(["x", "y"], 45)},
                        index=[f"M{i}" for i in range(90)])
    shift = meta["tissue"].map({"a": 0.0, "b": 5.0, "c": -3.0}) + meta["source"].map({"x": 0.0, "y": 2.0})
    values = pd.DataFrame({"g": shift + rng.normal(size=90), "h": rng.normal(size=90)}, index=meta.index)
    cleaned = clusters.remove_effects(values, meta, ["tissue", "source"])
    assert cleaned.groupby(meta["tissue"]).mean().abs().max().max() < 1e-9
    assert cleaned.groupby(meta["source"]).mean().abs().max().max() < 1e-9
    assert cleaned["g"].std() < values["g"].std()


def test_components_shape():
    table, _ = blobs()
    points = clusters.components(table, n_genes=30, n_components=5)
    assert points.shape == (80, 5)


def test_well_separated_groups_are_stable_and_found():
    table, truth = blobs()
    labels, agreement = clusters.stability(lambda rows: clusters.components(table.iloc[rows], n_components=5),
                                           len(table), k=2, n_subsamples=20)
    assert agreement.mean() > 0.95
    assert clusters.verdict(agreement.mean()) == "stable"
    assert (labels == truth).all() or (labels == 1 - truth).all()


def test_data_without_groups_is_less_stable_than_real_groups():
    rng = np.random.default_rng(1)
    cloud = rng.normal(size=(80, 10))                    # one round cloud: no groups, no preferred direction
    _, agreement = clusters.stability(lambda rows: cloud[rows], len(cloud), k=3, n_subsamples=20)
    assert agreement.mean() < 0.75


def test_reference_stability_is_reproducible_and_bounded():
    table, _ = blobs()
    points = clusters.components(table, n_components=5)
    first = clusters.reference_stability(points, k=2, n_reference=3, n_subsamples=5)
    second = clusters.reference_stability(points, k=2, n_reference=3, n_subsamples=5)
    np.testing.assert_array_equal(first, second)
    assert ((first > -0.5) & (first <= 1.0)).all()


def test_real_groups_beat_the_reference():
    table, _ = blobs()
    points = clusters.components(table, n_components=5)
    _, agreement = clusters.stability(lambda rows: clusters.components(table.iloc[rows], n_components=5),
                                      len(table), k=2, n_subsamples=20)
    reference = clusters.reference_stability(points, k=2, n_reference=5, n_subsamples=10)
    assert agreement.mean() > reference.max()


def test_verdict_thresholds():
    assert clusters.verdict(0.75) == "stable"
    assert clusters.verdict(0.60) == "in between"
    assert clusters.verdict(0.49) == "unstable"


def test_tissue_agreement_sees_groups_made_of_whole_tissues():
    """Three groups that each swallow two whole tissues: ARI is modest, but the groups follow tissue exactly."""
    from sklearn.metrics import adjusted_rand_score
    tissue = pd.Series(np.repeat(["a", "b", "c", "d", "e", "f"], 10))
    labels = np.repeat([0, 0, 1, 1, 2, 2], 10)
    share, expected = clusters.tissue_agreement(tissue, labels)
    assert share == pytest.approx(1.0)
    assert expected < 0.7
    assert adjusted_rand_score(tissue, labels) < 0.7


def test_tissue_agreement_ignores_tiny_tissues():
    tissue = pd.Series(["big"] * 10 + ["single"])
    labels = np.array([0] * 5 + [1] * 5 + [1])
    share, _ = clusters.tissue_agreement(tissue, labels)
    assert share == pytest.approx(0.5)                   # only the 10 'big' lines count; 5 are in its commonest group
