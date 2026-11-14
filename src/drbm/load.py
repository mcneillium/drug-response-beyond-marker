"""Read the raw files into tidy tables keyed by cell line (`model_id`).

Everything downstream uses four tables:

    response    one row per cell line: how it responded to one drug
    models      one row per cell line: tissue, cancer type
    pten        one row per cell line: is the PTEN gene altered?
    expression  cell lines x genes, log2(TPM + 1)
    protein     cell lines x proteins, log intensity (many gaps)

Nothing is filled or imputed here. A protein that was not detected in a cell
line stays missing.
"""

from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

GDSC_FILE = "GDSC2_fitted_dose_response_27Oct23.xlsx"
TRUNCATING = {"frameshift", "nonsense", "ess_splice", "stop_lost", "start_lost"}


def load_screen(raw: Path) -> pd.DataFrame:
    """The whole GDSC2 screen, cached as parquet because the Excel file is slow to read."""
    cache = Path(raw).parent / "gdsc2.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    screen = pd.read_excel(Path(raw) / GDSC_FILE)
    keep = ["SANGER_MODEL_ID", "CELL_LINE_NAME", "TCGA_DESC", "DRUG_ID", "DRUG_NAME",
            "PUTATIVE_TARGET", "PATHWAY_NAME", "LN_IC50", "AUC", "RMSE", "Z_SCORE"]
    screen = screen[keep].rename(columns={"SANGER_MODEL_ID": "model_id"})
    screen.to_parquet(cache, index=False)
    return screen


def load_response(raw: Path, drug: str) -> pd.DataFrame:
    """One row per cell line for one drug. AUC: 1 = no effect, lower = more sensitive.

    A few drugs were screened twice under different drug ids. The screen with
    more cell lines is used, so each cell line appears once.
    """
    screen = load_screen(raw)
    rows = screen[screen["DRUG_NAME"] == drug]
    if rows.empty:
        raise KeyError(f"Drug '{drug}' is not in the screen")
    largest = rows.groupby("DRUG_ID")["model_id"].nunique().idxmax()
    rows = rows[rows["DRUG_ID"] == largest].drop_duplicates("model_id")
    out = rows.set_index("model_id")[["CELL_LINE_NAME", "AUC", "LN_IC50", "RMSE"]]
    out.columns = ["cell_line", "auc", "ln_ic50", "curve_rmse"]
    return out.sort_index()


def load_models(raw: Path) -> pd.DataFrame:
    """Cell lines only (no organoids), with tissue and cancer type."""
    models = pd.read_csv(Path(raw) / "model_list_latest.csv.gz", low_memory=False)
    models = models[models["model_type"] == "Cell Line"]
    keep = ["model_id", "model_name", "tissue", "cancer_type", "msi_status",
            "mutation_data", "cnv_data", "growth_properties"]
    return models[keep].set_index("model_id").sort_index()


def load_pten(raw: Path, models: pd.DataFrame) -> pd.DataFrame:
    """PTEN status for every cell line that has both mutation and copy number data.

    altered     any cancer-driver mutation in PTEN, or both copies deleted
    truncating  a mutation that cuts the protein short, or both copies deleted
                (a stricter definition, used as a sensitivity check)

    A cell line without mutation or copy number data gets no status at all. It
    is not assumed to be normal.
    """
    raw = Path(raw)
    with zipfile.ZipFile(raw / "mutations_summary_20230202.zip") as archive:
        mutations = pd.read_csv(archive.open(archive.namelist()[0]))
    mutations = mutations[(mutations["gene_symbol"] == "PTEN") & (mutations["cancer_driver"] == True)]  # noqa: E712
    any_driver = set(mutations["model_id"])
    truncating = set(mutations.loc[mutations["effect"].isin(TRUNCATING), "model_id"])

    with zipfile.ZipFile(raw / "cnv_summary_20230303.zip") as archive:
        chunks = pd.read_csv(archive.open(archive.namelist()[0]), usecols=["model_id", "symbol", "cn_category"],
                             chunksize=1_000_000)
        copy_number = pd.concat(chunk[chunk["symbol"] == "PTEN"] for chunk in chunks)
    deleted = set(copy_number.loc[copy_number["cn_category"] == "Deletion", "model_id"])
    has_copy_number = set(copy_number["model_id"])

    known = models.index[(models["mutation_data"] == True) & models.index.isin(has_copy_number)]  # noqa: E712
    table = pd.DataFrame(index=known)
    table["pten_mutation"] = table.index.isin(any_driver)
    table["pten_deletion"] = table.index.isin(deleted)
    table["pten_altered"] = table["pten_mutation"] | table["pten_deletion"]
    table["pten_truncating"] = table.index.isin(truncating) | table["pten_deletion"]
    return table.sort_index()


def load_expression(raw: Path) -> tuple[pd.DataFrame, pd.Series]:
    """Gene expression as cell lines x genes, in log2(TPM + 1), and the sequencing source.

    Some cell lines were sequenced by both the Sanger and the Broad institutes.
    The Sanger run is used where there is one, the Broad run otherwise. The
    source is returned so it can be checked as a possible technical confounder.
    """
    raw = Path(raw)
    cache = raw.parent / "expression_log2tpm.parquet"
    source_cache = raw.parent / "expression_source.parquet"
    if cache.exists() and source_cache.exists():
        return pd.read_parquet(cache), pd.read_parquet(source_cache)["source"]

    with zipfile.ZipFile(raw / "rnaseq_merged_20250117.zip") as archive:
        with archive.open("rnaseq_merged_rsem_tpm_20250117.csv") as handle:
            header = pd.read_csv(handle, nrows=3, header=None)
        with archive.open("rnaseq_merged_rsem_tpm_20250117.csv") as handle:
            body = pd.read_csv(handle, skiprows=4, header=None, low_memory=False)
    model_ids = header.iloc[0, 3:].astype(str).to_numpy()
    sources = header.iloc[2, 3:].astype(str).to_numpy()
    symbols = body.iloc[:, 0].astype(str).to_numpy()
    # A small number of cells in the published file hold two numbers separated
    # by a space. They cannot be trusted, so they become missing; genes with
    # any missing value are dropped below.
    numeric = body.iloc[:, 3:].apply(pd.to_numeric, errors="coerce")
    values = numeric.to_numpy(dtype="float32")
    del body, numeric

    # One column per cell line: Sanger first, then Broad.
    order = np.argsort([0 if s == "Sanger" else 1 for s in sources], kind="stable")
    seen, keep = set(), []
    for position in order:
        if model_ids[position] not in seen:
            seen.add(model_ids[position])
            keep.append(position)
    keep = np.array(sorted(keep))
    matrix = pd.DataFrame(np.log2(values[:, keep] + 1.0).T, index=model_ids[keep], columns=symbols)
    matrix = matrix.loc[:, ~matrix.columns.duplicated()]
    matrix = matrix.loc[:, matrix.notna().all(axis=0)]
    matrix.index.name = "model_id"
    source = pd.Series(sources[keep], index=matrix.index, name="source")
    matrix.to_parquet(cache)
    source.to_frame().to_parquet(source_cache)
    return matrix, source


def load_protein(raw: Path) -> pd.DataFrame:
    """Protein abundance as cell lines x proteins (gene symbols). Undetected = NaN."""
    raw = Path(raw)
    cache = raw.parent / "protein.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    with zipfile.ZipFile(raw / "Proteomics_20221214.zip") as archive:
        with archive.open("Protein_matrix_averaged_20221214.tsv") as handle:
            table = pd.read_csv(handle, sep="\t", header=None, low_memory=False)
    symbols = table.iloc[1, 2:].astype(str).to_numpy()
    body = table.iloc[3:]
    matrix = pd.DataFrame(body.iloc[:, 2:].to_numpy(dtype="float32"), index=body.iloc[:, 1].astype(str).to_numpy(),
                          columns=symbols)
    matrix = matrix.loc[:, ~matrix.columns.duplicated()]
    matrix.index.name = "model_id"
    matrix.to_parquet(cache)
    return matrix
