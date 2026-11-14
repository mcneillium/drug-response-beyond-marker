import pandas as pd

from drbm import dataset


def test_group_rare_merges_small_groups_only():
    labels = pd.Series(["a"] * 5 + ["b"] * 2 + ["c"])
    grouped = dataset.group_rare(labels, minimum=3)
    assert grouped.tolist() == ["a"] * 5 + ["Other"] * 3
