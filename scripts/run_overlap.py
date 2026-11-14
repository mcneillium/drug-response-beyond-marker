"""Experiment 2: partial overlap between gene expression and protein data.

    python scripts/run_overlap.py

Protein data is kept for a share of cell lines and hidden for the rest, then
four strategies are compared (see src/drbm/overlap.py and docs/protocol.md).
The hiding is simulated: in the published data 97% of these cell lines have
protein data.

Outputs (results/)
    raw/overlap_predictions.parquet   held-out predictions, one row per cell line, repeat and share
    raw/overlap_fits.parquet          ridge strength, blend weight and training size of every fit
    table_overlap.csv / .md           R2 by share and strategy
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from drbm import PRIMARY_DRUG, dataset, folds, overlap, scores

ROOT = Path(__file__).resolve().parents[1]
SHARES = (0.97, 0.75, 0.50, 0.25)
LABELS = {"expression_only": "Expression only", "complete_cases": "Complete cases",
          "fill_gaps": "Fill the gaps", "late_fusion": "Late fusion"}


def run(data: dataset.Dataset, drug_folds: list) -> tuple[pd.DataFrame, pd.DataFrame]:
    predictions, fits = [], []
    availability = {}
    for repeat, fold, train, test in folds.iterate(drug_folds, data.ids):
        y_train = data.y[train]
        expr_train, expr_test = data.expression.loc[train], data.expression.loc[test]
        prot_train, prot_test = data.protein.loc[train], data.protein.loc[test]
        baseline, baseline_alpha = overlap.expression_only(y_train, expr_train, expr_test)
        expression_half = overlap.expression_part(y_train, expr_train, expr_test, seed=repeat)
        for share in SHARES:
            if (repeat, share) not in availability:
                availability[(repeat, share)] = overlap.keep_protein(data.has_protein, share, seed=repeat)
            avail = availability[(repeat, share)]
            avail_train, avail_test = avail[train], avail[test]
            complete, complete_alpha, n_complete = overlap.complete_cases(
                y_train, expr_train, expr_test, prot_train, prot_test, avail_train, avail_test)
            filled, filled_alpha = overlap.fill_gaps(
                y_train, expr_train, expr_test, prot_train, prot_test, avail_train, avail_test)
            fused, weight = overlap.late_fusion(
                y_train, expr_train, expr_test, prot_train, prot_test, avail_train, avail_test, seed=repeat,
                expression=expression_half)
            predictions.append(pd.DataFrame({
                "model_id": test, "repeat": repeat, "fold": fold, "share": share, "y": data.y[test].to_numpy(),
                "has_protein": avail_test.to_numpy(), "expression_only": baseline.to_numpy(),
                "complete_cases": complete.to_numpy(), "fill_gaps": filled.to_numpy(), "late_fusion": fused.to_numpy(),
            }))
            fits.append({"repeat": repeat, "fold": fold, "share": share, "n_train": len(train),
                         "n_train_with_protein": int(avail_train.sum()), "n_train_complete_cases": n_complete,
                         "alpha_expression_only": baseline_alpha, "alpha_complete_cases": complete_alpha,
                         "alpha_fill_gaps": filled_alpha, "late_fusion_weight_on_protein": weight})
        print(f"  repeat {repeat} fold {fold} done", flush=True)
    return pd.concat(predictions, ignore_index=True), pd.DataFrame(fits)


def summarise(predictions: pd.DataFrame, fits: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for share, at_share in predictions.groupby("share"):
        for strategy in overlap.STRATEGIES:
            per_repeat = []
            for _, block in at_share.groupby("repeat"):
                both, without = block[block["has_protein"]], block[~block["has_protein"]]
                covers_all = block[strategy].notna().all()
                per_repeat.append({
                    "r2_all": scores.r2(block["y"], block[strategy]) if covers_all else np.nan,
                    "r2_with_protein": scores.r2(both["y"], both[strategy]),
                    "r2_without_protein": scores.r2(without["y"], without[strategy])
                    if covers_all and len(without) >= 20 else np.nan,
                    "n_predicted": int(block[strategy].notna().sum()), "n_with_protein": len(both),
                })
            table = pd.DataFrame(per_repeat)
            share_fits = fits[fits["share"] == share]
            rows.append({
                "share_with_protein": share, "strategy": strategy, "label": LABELS[strategy],
                "cell_lines_predicted": int(table["n_predicted"].mean()),
                "training_lines_used": int(round(share_fits["n_train_complete_cases"].mean()
                                                 if strategy == "complete_cases" else share_fits["n_train"].mean())),
                **{f"{name}_mean": table[name].mean() for name in ("r2_all", "r2_with_protein", "r2_without_protein")},
                **{f"{name}_sd": table[name].std(ddof=1) for name in ("r2_all", "r2_with_protein", "r2_without_protein")},
                "mean_weight_on_protein": share_fits["late_fusion_weight_on_protein"].mean()
                if strategy == "late_fusion" else np.nan,
            })
    return pd.DataFrame(rows).sort_values(["share_with_protein", "strategy"], ascending=[False, True])


def to_markdown(table: pd.DataFrame, n_lines: int) -> str:
    def cell(mean: float, sd: float) -> str:
        return "n/a" if np.isnan(mean) else f"{mean:.3f} ({sd:.3f})"
    lines = [f"{PRIMARY_DRUG}, {n_lines} cell lines. Held-out R², mean over 5 repeats (standard deviation).",
             "Protein data hidden at random for part of the cell lines (simulated).", "",
             "| Share with protein | Strategy | Training lines used | Lines predicted | R², all lines | "
             "R², lines with protein | R², lines without protein |",
             "| ---: | --- | ---: | ---: | ---: | ---: | ---: |"]
    for _, row in table.iterrows():
        lines.append(f"| {row['share_with_protein']:.2f} | {row['label']} | {row['training_lines_used']} | "
                     f"{row['cell_lines_predicted']} | {cell(row['r2_all_mean'], row['r2_all_sd'])} | "
                     f"{cell(row['r2_with_protein_mean'], row['r2_with_protein_sd'])} | "
                     f"{cell(row['r2_without_protein_mean'], row['r2_without_protein_sd'])} |")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "raw")
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    raw_dir = args.results / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    data = dataset.assemble(args.data, PRIMARY_DRUG)
    drug_folds = folds.load_folds(ROOT / "splits" / "folds_v1.json")[PRIMARY_DRUG]
    predictions_path, fits_path = raw_dir / "overlap_predictions.parquet", raw_dir / "overlap_fits.parquet"
    if not predictions_path.exists():
        predictions, fits = run(data, drug_folds)
        fits.to_parquet(fits_path, index=False)
        predictions.to_parquet(predictions_path, index=False)
    predictions, fits = pd.read_parquet(predictions_path), pd.read_parquet(fits_path)
    table = summarise(predictions, fits)
    table.to_csv(args.results / "table_overlap.csv", index=False)
    (args.results / "table_overlap.md").write_text(to_markdown(table, len(data.ids)))
    print(to_markdown(table, len(data.ids)))


if __name__ == "__main__":
    main()
