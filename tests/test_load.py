"""The loaders, checked on tiny hand-made files with the same layout as the real ones."""

import zipfile

import numpy as np
import pandas as pd
import pytest

from drbm import load


@pytest.fixture
def raw(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    screen = pd.DataFrame({
        "model_id": ["M1", "M2", "M3", "M1", "M1"],
        "CELL_LINE_NAME": ["one", "two", "three", "one", "one"],
        "TCGA_DESC": "X", "DRUG_ID": [10, 10, 10, 11, 20], "DRUG_NAME": ["DrugA", "DrugA", "DrugA", "DrugA", "DrugB"],
        "PUTATIVE_TARGET": "T", "PATHWAY_NAME": "P", "LN_IC50": [1.0, 2.0, 3.0, 9.0, 5.0],
        "AUC": [0.5, 0.7, 0.9, 0.1, 0.6], "RMSE": 0.05, "Z_SCORE": 0.0,
    })
    screen.to_parquet(tmp_path / "gdsc2.parquet", index=False)

    mutations = pd.DataFrame({
        "model_id": ["M1", "M2", "M3", "M4"],
        "gene_symbol": ["PTEN", "PTEN", "TP53", "PTEN"],
        "cancer_driver": [True, True, True, False],
        "effect": ["frameshift", "missense", "nonsense", "nonsense"],
    })
    with zipfile.ZipFile(raw / "mutations_summary_20230202.zip", "w") as archive:
        archive.writestr("mutations_summary.csv", mutations.to_csv(index=False))
    copy_number = pd.DataFrame({
        "model_id": ["M1", "M2", "M3", "M4", "M5", "M5"],
        "symbol": ["PTEN", "PTEN", "PTEN", "PTEN", "PTEN", "TP53"],
        "cn_category": ["Neutral", "Neutral", "Neutral", "Neutral", "Deletion", "Neutral"],
    })
    with zipfile.ZipFile(raw / "cnv_summary_20230303.zip", "w") as archive:
        archive.writestr("cnv_summary.csv", copy_number.to_csv(index=False))
    models = pd.DataFrame({"mutation_data": [True, True, True, True, True, True, False]},
                          index=pd.Index(["M1", "M2", "M3", "M4", "M5", "M6", "M7"], name="model_id"))
    return raw, models


def test_response_uses_the_larger_screen_and_one_row_per_cell_line(raw):
    raw_dir, _ = raw
    response = load.load_response(raw_dir, "DrugA")
    assert response.index.tolist() == ["M1", "M2", "M3"]
    assert response.loc["M1", "auc"] == 0.5             # from screen 10 (3 lines), not screen 11 (1 line)
    with pytest.raises(KeyError):
        load.load_response(raw_dir, "Not a drug")


def test_pten_definitions(raw):
    raw_dir, models = raw
    pten = load.load_pten(raw_dir, models)
    # M1: truncating driver mutation.  M2: driver mutation, not truncating.
    # M3: mutation is in another gene.  M4: PTEN mutation that is not a driver.
    # M5: both copies deleted.
    assert pten.loc["M1", ["pten_altered", "pten_truncating"]].tolist() == [True, True]
    assert pten.loc["M2", ["pten_altered", "pten_truncating"]].tolist() == [True, False]
    assert pten.loc["M3", ["pten_altered", "pten_truncating"]].tolist() == [False, False]
    assert pten.loc["M4", ["pten_altered", "pten_truncating"]].tolist() == [False, False]
    assert pten.loc["M5", ["pten_altered", "pten_truncating"]].tolist() == [True, True]


def test_pten_status_is_not_assumed_when_data_is_missing(raw):
    raw_dir, models = raw
    pten = load.load_pten(raw_dir, models)
    assert "M6" not in pten.index                        # no copy number data
    assert "M7" not in pten.index                        # no mutation data
