"""Experiment 6 (exploratory, added after seeing results): the pathway genes, forced in.

    python scripts/run_pathway_check.py

The 2,000-most-variable-genes filter removed 10 of the 12 declared pathway
genes, including PTEN itself, so the main models never saw them. This script
asks what that cost, on the same folds as everything else:

    P1  the expression level of the PTEN gene alone (one straight line)
    P2  PTEN mutation status + PTEN expression level
    P3  PTEN status + background + the 12 pathway genes (ridge regression)

Compare P1 and P2 with L1 (PTEN status), and P3 with L3 (PTEN status +
background). Read docs/protocol.md, amendment 8, before quoting these numbers.

Output: results/table_pathway_check.csv / .md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from drbm import PRIMARY_DRUG, dataset, folds, models, scores

ROOT = Path(__file__).resolve().parents[1]
LABELS = {"L1": "PTEN status (from Experiment 1)", "P1": "PTEN gene expression level alone",
          "P2": "PTEN status + PTEN expression level", "L3": "PTEN status + background (from Experiment 1)",
          "P3": "PTEN status + background + 12 pathway genes"}


def fit(y_train: pd.Series, train: pd.DataFrame, test: pd.DataFrame) -> pd.Series:
    train_x, test_x = models.standardise(train, test)
    prediction, _ = models.ridge(train_x, y_train.to_numpy(), test_x)
    return pd.Series(prediction, index=test.index)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "raw")
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    args = parser.parse_args()

    data = dataset.assemble(args.data, PRIMARY_DRUG)
    drug_folds = folds.load_folds(ROOT / "splits" / "folds_v1.json")[PRIMARY_DRUG]
    genes = [g for g in models.PATHWAY_GENES if g in data.expression.columns]
    missing = sorted(set(models.PATHWAY_GENES) - set(genes))
    pathway = data.expression[genes].astype(float)

    parts = []
    for repeat, fold, train, test in folds.iterate(drug_folds, data.ids):
        y_train, meta_train, meta_test = data.y[train], data.meta.loc[train], data.meta.loc[test]
        status_train, status_test = models.pten_column(meta_train), models.pten_column(meta_test)
        bg_train, bg_test = models.background(meta_train, meta_test)
        block = pd.DataFrame({"repeat": repeat, "y": data.y[test]}, index=test)
        block["P1"] = fit(y_train, pathway.loc[train, ["PTEN"]], pathway.loc[test, ["PTEN"]])
        block["P2"] = fit(y_train, pd.concat([status_train, pathway.loc[train, ["PTEN"]]], axis=1),
                          pd.concat([status_test, pathway.loc[test, ["PTEN"]]], axis=1))
        block["P3"] = fit(y_train, pd.concat([status_train, bg_train, pathway.loc[train]], axis=1),
                          pd.concat([status_test, bg_test, pathway.loc[test]], axis=1))
        parts.append(block.rename_axis("model_id").reset_index())
    new = pd.concat(parts, ignore_index=True)
    old = pd.read_parquet(args.results / "raw" / "oof_ipatasertib.parquet")[["model_id", "repeat", "L1", "L3"]]
    oof = new.merge(old, on=["model_id", "repeat"], validate="one_to_one")

    rows = []
    for step in ("L1", "P1", "P2", "L3", "P3"):
        per_repeat = pd.DataFrame([{
            "r2": scores.r2(b["y"].to_numpy(), b[step].to_numpy()),
            "spearman": scores.spearman(b["y"].to_numpy(), b[step].to_numpy()),
            "within_type_spearman": scores.within_group_spearman(
                b.set_index("model_id")["y"], b.set_index("model_id")[step], data.meta["cancer_type"]),
        } for _, b in oof.groupby("repeat")])
        rows.append({"step": step, "label": LABELS[step], **{f"{c}_mean": per_repeat[c].mean() for c in per_repeat},
                     **{f"{c}_sd": per_repeat[c].std(ddof=1) for c in per_repeat}})
    table = pd.DataFrame(rows)

    wide = oof.pivot(index="model_id", columns="repeat")
    observed = wide["y"].iloc[:, 0].to_numpy()
    gains = []
    for better, simpler in (("P1", "L1"), ("P2", "L1"), ("P3", "L3")):
        point, low, high = scores.repeated_gain(observed, wide[better].to_numpy(), wide[simpler].to_numpy())
        gains.append({"better": better, "simpler": simpler, "gain_r2": point, "ci_low": low, "ci_high": high})
    gains = pd.DataFrame(gains)
    table.to_csv(args.results / "table_pathway_check.csv", index=False)
    gains.to_csv(args.results / "table_pathway_check_gains.csv", index=False)

    lines = ["EXPLORATORY: added after seeing the results of Experiments 1 to 5.",
             f"{PRIMARY_DRUG}, {len(data.ids)} cell lines, same folds. Mean over 5 repeats (standard deviation).",
             f"Pathway genes found in the expression data: {len(genes)} of {len(models.PATHWAY_GENES)}"
             + (f" (missing: {', '.join(missing)})." if missing else "."), "",
             "| Step | Inputs | R² | Spearman | Within-cancer-type Spearman |", "| --- | --- | ---: | ---: | ---: |"]
    for _, row in table.iterrows():
        lines.append(f"| {row['step']} | {row['label']} | {row['r2_mean']:.3f} ({row['r2_sd']:.3f}) | "
                     f"{row['spearman_mean']:.3f} ({row['spearman_sd']:.3f}) | "
                     f"{row['within_type_spearman_mean']:.3f} ({row['within_type_spearman_sd']:.3f}) |")
    lines += ["", "| Better | Simpler | Gain in R² | 95% interval |", "| --- | --- | ---: | --- |"]
    for _, row in gains.iterrows():
        lines.append(f"| {row['better']} | {row['simpler']} | {row['gain_r2']:+.3f} | "
                     f"{row['ci_low']:+.3f} to {row['ci_high']:+.3f} |")
    (args.results / "table_pathway_check.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
