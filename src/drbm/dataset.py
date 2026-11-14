"""Assemble one drug's analysis tables: response, cell line details, expression, protein.

Only cell lines with all three of a drug response, a known PTEN status and
gene expression are kept. Protein data is attached where it exists; cell lines
without it keep a row of missing values.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from drbm import load

RARE_CANCER_TYPE = 10      # cancer types with fewer cell lines are grouped as "Other"
RARE_TISSUE = 5            # tissues with fewer cell lines are grouped for fold stratification


@dataclass
class Dataset:
    drug: str
    y: pd.Series               # AUC per cell line; lower = more sensitive
    meta: pd.DataFrame         # tissue, cancer_type, cancer_group, source, pten_altered, pten_truncating
    expression: pd.DataFrame   # cell lines x genes, log2(TPM + 1)
    protein: pd.DataFrame      # cell lines x proteins; all-missing rows where there is no protein data

    @property
    def ids(self) -> pd.Index:
        return self.y.index

    @property
    def has_protein(self) -> pd.Series:
        return self.protein.notna().any(axis=1)


def group_rare(labels: pd.Series, minimum: int, other: str = "Other") -> pd.Series:
    """Replace labels that occur fewer than `minimum` times with one shared label."""
    counts = labels.value_counts()
    return labels.where(labels.map(counts) >= minimum, other)


def assemble(raw: Path, drug: str) -> Dataset:
    raw = Path(raw)
    models = load.load_models(raw)
    pten = load.load_pten(raw, models)
    response = load.load_response(raw, drug)
    expression, source = load.load_expression(raw)
    protein = load.load_protein(raw)

    ids = response.index.intersection(pten.index).intersection(expression.index).sort_values()
    meta = models.loc[ids, ["model_name", "tissue", "cancer_type"]].copy()
    meta["cancer_group"] = group_rare(meta["cancer_type"], RARE_CANCER_TYPE)
    meta["tissue_group"] = group_rare(meta["tissue"], RARE_TISSUE)
    meta["source"] = source.loc[ids]
    meta["pten_altered"] = pten.loc[ids, "pten_altered"].astype(bool)
    meta["pten_truncating"] = pten.loc[ids, "pten_truncating"].astype(bool)
    return Dataset(
        drug=drug,
        y=response.loc[ids, "auc"].astype(float).rename("auc"),
        meta=meta,
        expression=expression.loc[ids],
        protein=protein.reindex(ids),
    )
