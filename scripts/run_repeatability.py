"""Experiment 7 (exploratory, added after seeing results): how repeatable is the response?

    python scripts/run_repeatability.py

Finds every drug that was measured twice on the same cell lines, either under
two drug ids within GDSC2 or in both GDSC1 and GDSC2, and reports how well the
two sets of AUC values agree. Read docs/protocol.md, amendment 10, first.

Outputs (results/)
    table_repeatability.csv       every repeated drug
    table_repeatability.md        the two AKT inhibitors, and a summary of the rest
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from drbm import repeats

ROOT = Path(__file__).resolve().parents[1]
FILES = {"GDSC1": "GDSC1_fitted_dose_response_27Oct23.xlsx", "GDSC2": "GDSC2_fitted_dose_response_27Oct23.xlsx"}
COLUMNS = ["SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET", "MIN_CONC", "MAX_CONC", "AUC"]
HIGHLIGHT = {"Uprosertib": "within GDSC2", "MK-2206": "GDSC1 against GDSC2"}


def load_screens(raw: Path) -> pd.DataFrame:
    cache = raw.parent / "screens_for_repeatability.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    parts = []
    for screen, name in FILES.items():
        table = pd.read_excel(raw / name, usecols=COLUMNS).rename(columns={"SANGER_MODEL_ID": "model_id"})
        parts.append(table.assign(screen=screen))
    screens = pd.concat(parts, ignore_index=True)
    screens.to_parquet(cache, index=False)
    return screens


def one_measurement(screens: pd.DataFrame, screen: str, drug_id: int) -> tuple[pd.Series, str]:
    rows = screens[(screens["screen"] == screen) & (screens["DRUG_ID"] == drug_id)].drop_duplicates("model_id")
    doses = f"{rows['MIN_CONC'].median():g} to {rows['MAX_CONC'].median():g}"
    return rows.set_index("model_id")["AUC"], doses


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "raw")
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    screens = load_screens(args.data)

    pairs = []
    # Within GDSC2: one drug name under two drug ids. Each id is a separate screening run.
    second = screens[screens["screen"] == "GDSC2"]
    for drug, rows in second.groupby("DRUG_NAME"):
        ids = rows.groupby("DRUG_ID")["model_id"].nunique().sort_values(ascending=False).index.tolist()
        if len(ids) >= 2:
            pairs.append((drug, "within GDSC2", ("GDSC2", ids[0]), ("GDSC2", ids[1])))
    # Across screens: the same drug in GDSC1 and GDSC2 (largest run in each).
    first = screens[screens["screen"] == "GDSC1"]
    for drug in sorted(set(first["DRUG_NAME"]) & set(second["DRUG_NAME"])):
        largest = [rows[rows["DRUG_NAME"] == drug].groupby("DRUG_ID")["model_id"].nunique().idxmax()
                   for rows in (first, second)]
        pairs.append((drug, "GDSC1 against GDSC2", ("GDSC1", largest[0]), ("GDSC2", largest[1])))

    rows = []
    for drug, kind, a, b in pairs:
        first_auc, first_doses = one_measurement(screens, *a)
        second_auc, second_doses = one_measurement(screens, *b)
        highlighted = HIGHLIGHT.get(drug) == kind
        result = repeats.pair_agreement(first_auc, second_auc, n_resamples=2000 if highlighted else 0)
        rows.append({"drug": drug, "kind": kind, "first": f"{a[0]} id {a[1]}", "second": f"{b[0]} id {b[1]}",
                     "doses_first_uM": first_doses, "doses_second_uM": second_doses,
                     "same_dose_range": first_doses == second_doses, "highlighted": highlighted, **result})
    table = pd.DataFrame(rows).dropna(subset=["pearson"]).sort_values(["kind", "drug"])
    table.to_csv(args.results / "table_repeatability.csv", index=False)

    replication = pd.read_csv(args.results / "table_replication.csv").set_index("drug")
    lines = ["EXPLORATORY: added after seeing the results of Experiments 1 to 6.",
             "Agreement between two measurements of the same drug on the same cell lines (AUC).", "",
             "| Drug | Repeat | Cell lines | Dose ranges (µM) | Pearson (95% interval) | Spearman | "
             "Best model R² here | Share of the ceiling |", "| --- | --- | ---: | --- | --- | ---: | ---: | ---: |"]
    for _, row in table[table["highlighted"]].iterrows():
        model = float(replication.loc[row["drug"], "r2_L4_expression"])
        lines.append(f"| {row['drug']} | {row['kind']} | {row['n_shared']} | {row['doses_first_uM']} and "
                     f"{row['doses_second_uM']} | {row['pearson']:.2f} ({row['pearson_low']:.2f} to "
                     f"{row['pearson_high']:.2f}) | {row['spearman']:.2f} | {model:.3f} | "
                     f"{repeats.share_of_ceiling(model, row['pearson']):.0%} |")
    lines += ["", "\"Best model R² here\" is the elastic net with expression (L4) for that drug, from Experiment 5.",
              "\"Share of the ceiling\" is that R² divided by the Pearson correlation.", "",
              "| All repeated drugs | Number | Median Pearson | Lowest | Highest |", "| --- | ---: | ---: | ---: | ---: |"]
    for kind, block in table.groupby("kind"):
        lines.append(f"| {kind} | {len(block)} | {block['pearson'].median():.2f} | {block['pearson'].min():.2f} | "
                     f"{block['pearson'].max():.2f} |")
        same = block[block["same_dose_range"]]
        if 0 < len(same) < len(block):
            lines.append(f"| {kind}, same dose range only | {len(same)} | {same['pearson'].median():.2f} | "
                         f"{same['pearson'].min():.2f} | {same['pearson'].max():.2f} |")
    (args.results / "table_repeatability.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
