"""Small synthetic tables shared by the tests. No real data is needed to run them."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def toy():
    """120 'cell lines', 60 'genes', 30 'proteins'. Gene g0 drives the response."""
    rng = np.random.default_rng(0)
    ids = pd.Index([f"M{i:03d}" for i in range(120)], name="model_id")
    expression = pd.DataFrame(rng.normal(size=(120, 60)), index=ids, columns=[f"g{j}" for j in range(60)])
    protein = pd.DataFrame(rng.normal(size=(120, 30)), index=ids, columns=[f"p{j}" for j in range(30)])
    y = pd.Series(0.5 + 0.1 * expression["g0"] + 0.02 * rng.normal(size=120), index=ids, name="auc")
    meta = pd.DataFrame({
        "tissue": np.where(np.arange(120) % 3 == 0, "Lung", "Skin"),
        "cancer_group": np.where(np.arange(120) % 3 == 0, "Lung cancer", "Melanoma"),
        "source": np.where(np.arange(120) % 2 == 0, "Sanger", "Broad"),
        "pten_altered": np.arange(120) % 5 == 0,
        "pten_truncating": np.arange(120) % 10 == 0,
    }, index=ids)
    meta["cancer_type"] = meta["cancer_group"]
    meta["tissue_group"] = meta["tissue"]
    return {"ids": ids, "expression": expression, "protein": protein, "y": y, "meta": meta}
