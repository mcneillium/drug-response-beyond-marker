"""Experiments 1, 1b, 4 and 5: the baseline ladder, the leak, explanation stability, replication.

    python scripts/run_ladder.py

Follows docs/protocol.md. Held-out predictions are saved per drug under
results/raw/, so a finished drug is not refitted on a re-run.

Outputs (results/)
    table_ladder.csv / .md          R2 and rank correlations for each ladder step
    table_ladder_gains.csv          gains over simpler steps, with 95% intervals
    table_leak.csv                  feature selection done wrongly and correctly
    table_feature_stability.csv     how often each gene is selected / in the top 20
    table_explanation_summary.csv   overlap of top-20 lists; where pathway genes rank
    table_replication.csv           the ladder on the other AKT inhibitors
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from drbm import PRIMARY_DRUG, REPLICATION_DRUGS, dataset, folds, models, scores

ROOT = Path(__file__).resolve().parents[1]
STEPS = ["L0", "L1", "L2", "L3", "L4", "L5"]
STEP_LABELS = {
    "L0": "Training mean", "L1": "PTEN status", "L2": "Background factors",
    "L3": "PTEN + background", "L4": "+ expression (elastic net)", "L5": "+ expression (boosting)",
    "L1_strict": "PTEN status, strict definition", "L3_strict": "PTEN strict + background",
}
TOP_K = 20
LEAK_GENES = 200


def run_drug(data: dataset.Dataset, drug_folds: list, full: bool) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Held-out predictions for every ladder step, plus per-fit weights and SHAP importance.

    `full` adds the boosting step and the strict-PTEN sensitivity check
    (primary drug only).
    """
    predictions, weights, importances = [], [], []
    for repeat, fold, train, test in folds.iterate(drug_folds, data.ids):
        y_train, meta_train, meta_test = data.y[train], data.meta.loc[train], data.meta.loc[test]
        expr_train, expr_test = data.expression.loc[train], data.expression.loc[test]
        block = pd.DataFrame({"repeat": repeat, "fold": fold, "y": data.y[test]}, index=test)
        block["L0"] = models.l0_mean(y_train, test)
        block["L1"] = models.l1_pten(y_train, meta_train, meta_test)
        block["L3"] = models.l3_pten_background(y_train, meta_train, meta_test)
        block["L4"], weight = models.l4_elastic_net(y_train, expr_train, expr_test, meta_train, meta_test, seed=repeat)
        weights.append(weight.rename("weight").rename_axis("feature").reset_index().assign(repeat=repeat, fold=fold))
        if full:
            block["L2"] = models.l2_background(y_train, meta_train, meta_test)
            block["L1_strict"] = models.l1_pten(y_train, meta_train, meta_test, "pten_truncating")
            block["L3_strict"] = models.l3_pten_background(y_train, meta_train, meta_test, "pten_truncating")
            block["L5"], importance = models.l5_boosting(y_train, expr_train, expr_test, meta_train, meta_test, seed=repeat)
            importances.append(
                importance.rename("mean_abs_shap").rename_axis("feature").reset_index().assign(repeat=repeat, fold=fold)
            )
        predictions.append(block.rename_axis("model_id").reset_index())
        print(f"  {data.drug}: repeat {repeat} fold {fold} done", flush=True)
    return (pd.concat(predictions, ignore_index=True), pd.concat(weights, ignore_index=True),
            pd.concat(importances, ignore_index=True) if importances else pd.DataFrame())


def score_steps(oof: pd.DataFrame, meta: pd.DataFrame, steps: list[str]) -> pd.DataFrame:
    """Mean and standard deviation over repeats of each score, per ladder step."""
    rows = []
    for step in steps:
        per_repeat = []
        for _, block in oof.groupby("repeat"):
            block = block.set_index("model_id")
            per_repeat.append({
                "r2": scores.r2(block["y"].to_numpy(), block[step].to_numpy()),
                "spearman": scores.spearman(block["y"].to_numpy(), block[step].to_numpy()),
                "within_type_spearman": scores.within_group_spearman(
                    block["y"], block[step], meta.loc[block.index, "cancer_type"]),
            })
        table = pd.DataFrame(per_repeat)
        rows.append({"step": step, "label": STEP_LABELS[step],
                     **{f"{name}_mean": table[name].mean() for name in table},
                     **{f"{name}_sd": table[name].std(ddof=1) for name in table}})
    return pd.DataFrame(rows)


def gains(oof: pd.DataFrame, pairs: list[tuple[str, str]]) -> pd.DataFrame:
    """Gain in R2 of one step over another: computed within each repeat, then averaged."""
    wide = oof.pivot(index="model_id", columns="repeat")
    observed = wide["y"].iloc[:, 0].to_numpy()
    rows = []
    for better, simpler in pairs:
        point, low, high = scores.repeated_gain(observed, wide[better].to_numpy(), wide[simpler].to_numpy())
        rows.append({"better": better, "simpler": simpler, "gain_r2": point, "ci_low": low, "ci_high": high,
                     "interval_excludes_zero": bool(low > 0 or high < 0)})
    return pd.DataFrame(rows)


def leak_experiment(data: dataset.Dataset, drug_folds: list) -> pd.DataFrame:
    """Experiment 1b: choose genes with all the data (wrong) or inside each fold (right)."""
    rng = np.random.default_rng(0)
    shuffled = pd.Series(rng.permutation(data.y.to_numpy()), index=data.ids)
    rows = []
    for label, y in (("real response", data.y), ("shuffled response", shuffled)):
        chosen_with_all_data = models.top_correlated(data.expression, y, LEAK_GENES)
        collected = {"wrong: genes chosen using all cell lines": [], "right: genes chosen inside each training fold": []}
        for repeat, fold, train, test in folds.iterate(drug_folds, data.ids):
            expr_train, expr_test = data.expression.loc[train], data.expression.loc[test]
            wrong = models.ridge_on_genes(y[train], expr_train, expr_test, chosen_with_all_data)
            inside = models.top_correlated(expr_train, y[train], LEAK_GENES)
            right = models.ridge_on_genes(y[train], expr_train, expr_test, inside)
            for name, prediction in zip(collected, (wrong, right)):
                collected[name].append(pd.DataFrame({"repeat": repeat, "y": y[test], "prediction": prediction}))
        for name, parts in collected.items():
            table = pd.concat(parts)
            per_repeat = [scores.r2(b["y"].to_numpy(), b["prediction"].to_numpy()) for _, b in table.groupby("repeat")]
            ranks = [scores.spearman(b["y"].to_numpy(), b["prediction"].to_numpy()) for _, b in table.groupby("repeat")]
            rows.append({"response": label, "method": name, "r2_mean": float(np.mean(per_repeat)),
                         "r2_sd": float(np.std(per_repeat, ddof=1)), "spearman_mean": float(np.mean(ranks)),
                         "spearman_sd": float(np.std(ranks, ddof=1))})
        print(f"  leak experiment: {label} done", flush=True)
    return pd.DataFrame(rows)


def feature_stability(weights: pd.DataFrame, importances: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Experiment 4: how repeatable are the selected genes and the top-20 lists?"""
    n_fits = weights.groupby(["repeat", "fold"]).ngroups
    selected = weights.groupby("feature").agg(times_selected=("weight", "size"), mean_weight=("weight", "mean"))
    selected["share_of_fits_selected"] = selected["times_selected"] / n_fits

    # A gene is only an input in the fits where it was among the 2,000 most
    # variable genes of that training fold. Importance is averaged over ALL
    # fits, counting zero where the gene was not an input, so a gene that
    # appears in three fits cannot outrank one that matters in all 25.
    importances = importances.copy()
    importances["rank"] = importances.groupby(["repeat", "fold"])["mean_abs_shap"].rank(ascending=False, method="first")
    top = importances[importances["rank"] <= TOP_K]
    in_top = top.groupby("feature").size().rename("times_in_top20")
    overall = importances.groupby("feature").agg(times_an_input=("mean_abs_shap", "size"),
                                                 total_shap=("mean_abs_shap", "sum"))
    overall["mean_abs_shap"] = overall.pop("total_shap") / n_fits
    table = overall.join(in_top).join(selected[["times_selected", "share_of_fits_selected", "mean_weight"]])
    table["times_in_top20"] = table["times_in_top20"].fillna(0).astype(int)
    table["share_of_fits_in_top20"] = table["times_in_top20"] / n_fits
    table["times_selected"] = table["times_selected"].fillna(0).astype(int)
    table["share_of_fits_selected"] = table["share_of_fits_selected"].fillna(0.0)
    table["overall_rank_by_shap"] = table["mean_abs_shap"].rank(ascending=False, method="first").astype(int)
    table = table.sort_values("overall_rank_by_shap").reset_index()

    top_sets = [set(block["feature"]) for _, block in top.groupby(["repeat", "fold"])]
    jaccard = [len(a & b) / len(a | b) for a, b in itertools.combinations(top_sets, 2)]
    selected_sets = [set(block["feature"]) for _, block in weights.groupby(["repeat", "fold"])]
    selected_jaccard = [len(a & b) / len(a | b) for a, b in itertools.combinations(selected_sets, 2)]
    rows = [
        {"quantity": "number of fits", "value": n_fits},
        {"quantity": "boosting: mean Jaccard overlap of top-20 lists", "value": float(np.mean(jaccard))},
        {"quantity": "boosting: features in the top 20 of at least 80% of fits", "value": int((table["share_of_fits_in_top20"] >= 0.8).sum())},
        {"quantity": "boosting: features that reach any top-20 list", "value": int((table["times_in_top20"] > 0).sum())},
        {"quantity": "features that are an input in every fit", "value": int((table["times_an_input"] == n_fits).sum())},
        {"quantity": "features that are an input in at least one fit", "value": len(table)},
        {"quantity": "elastic net: mean number of features selected per fit", "value": float(np.mean([len(s) for s in selected_sets]))},
        {"quantity": "elastic net: mean Jaccard overlap of selected sets", "value": float(np.mean(selected_jaccard))},
        {"quantity": "elastic net: features selected in at least 80% of fits", "value": int((table["share_of_fits_selected"] >= 0.8).sum())},
    ]
    for name in ["PTEN_status", *models.PATHWAY_GENES]:
        hit = table[table["feature"] == name]
        if hit.empty:
            rows.append({"quantity": f"sanity check: {name}",
                         "value": "never an input: not among the 2,000 most variable genes in any fit"})
        else:
            row = hit.iloc[0]
            rows.append({"quantity": f"sanity check: {name}",
                         "value": f"an input in {int(row['times_an_input'])}/{n_fits} fits; SHAP rank "
                                  f"{int(row['overall_rank_by_shap'])} of {len(table)}; in top 20 of "
                                  f"{int(row['times_in_top20'])}/{n_fits} fits; selected by elastic net in "
                                  f"{int(row['times_selected'])}/{n_fits}"})
    return table, pd.DataFrame(rows)


def ladder_markdown(table: pd.DataFrame, gain_table: pd.DataFrame, n_lines: int, drug: str) -> str:
    lines = [f"{drug}, {n_lines} cell lines. Held-out predictions, 5-fold cross-validation repeated 5 times.",
             "Mean over repeats (standard deviation in brackets).", "",
             "| Step | Inputs | R² | Spearman | Within-cancer-type Spearman |", "| --- | --- | ---: | ---: | ---: |"]
    for _, row in table.iterrows():
        lines.append(f"| {row['step']} | {row['label']} | {row['r2_mean']:.3f} ({row['r2_sd']:.3f}) | "
                     f"{row['spearman_mean']:.3f} ({row['spearman_sd']:.3f}) | "
                     f"{row['within_type_spearman_mean']:.3f} ({row['within_type_spearman_sd']:.3f}) |")
    lines += ["", "**Gains in R² (95% bootstrap interval over cell lines)**", "",
              "| Better | Simpler | Gain | 95% interval |", "| --- | --- | ---: | --- |"]
    for _, row in gain_table.iterrows():
        lines.append(f"| {row['better']} | {row['simpler']} | {row['gain_r2']:+.3f} | "
                     f"{row['ci_low']:+.3f} to {row['ci_high']:+.3f} |")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "raw")
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    raw_dir = args.results / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    drugs = [PRIMARY_DRUG, *REPLICATION_DRUGS]
    datasets = {drug: dataset.assemble(args.data, drug) for drug in drugs}
    folds_path = ROOT / "splits" / "folds_v1.json"
    if not folds_path.exists():
        folds.save_folds({drug: folds.make_folds(d.ids, d.meta["tissue_group"]) for drug, d in datasets.items()}, folds_path)
    all_folds = folds.load_folds(folds_path)

    outputs = {}
    for drug in drugs:
        data, full = datasets[drug], drug == PRIMARY_DRUG
        slug = drug.lower().replace("-", "")
        paths = {name: raw_dir / f"{name}_{slug}.parquet" for name in ("oof", "weights", "importances")}
        if not paths["oof"].exists():
            oof, weights, importances = run_drug(data, all_folds[drug], full)
            weights.to_parquet(paths["weights"], index=False)
            if full:
                importances.to_parquet(paths["importances"], index=False)
            oof.to_parquet(paths["oof"], index=False)
        outputs[drug] = {name: pd.read_parquet(path) for name, path in paths.items() if path.exists()}

    # ----- Experiment 1: the ladder, primary drug -----
    primary, data = outputs[PRIMARY_DRUG], datasets[PRIMARY_DRUG]
    ladder = score_steps(primary["oof"], data.meta, [*STEPS, "L1_strict", "L3_strict"])
    gain_table = gains(primary["oof"], [("L1", "L0"), ("L3", "L2"), ("L3", "L1"), ("L4", "L3"), ("L5", "L3"),
                                        ("L5", "L4"), ("L4", "L1"), ("L5", "L1")])
    ladder.to_csv(args.results / "table_ladder.csv", index=False)
    gain_table.to_csv(args.results / "table_ladder_gains.csv", index=False)
    (args.results / "table_ladder.md").write_text(ladder_markdown(ladder, gain_table, len(data.ids), PRIMARY_DRUG))

    # ----- Experiment 1b: the leak -----
    leak_path = args.results / "table_leak.csv"
    if not leak_path.exists():
        leak_experiment(data, all_folds[PRIMARY_DRUG]).to_csv(leak_path, index=False)

    # ----- Experiment 4: stability of explanations -----
    stability, summary = feature_stability(primary["weights"], primary["importances"])
    stability.to_csv(args.results / "table_feature_stability.csv", index=False)
    summary.to_csv(args.results / "table_explanation_summary.csv", index=False)

    # ----- Experiment 5: replication on the other AKT inhibitors -----
    n_fits = primary["weights"].groupby(["repeat", "fold"]).ngroups
    stable_primary = set(stability.loc[stability["share_of_fits_selected"] >= 0.8, "feature"])
    rows = []
    for drug in drugs:
        out, d = outputs[drug], datasets[drug]
        table = score_steps(out["oof"], d.meta, ["L1", "L3", "L4"]).set_index("step")
        gain = gains(out["oof"], [("L4", "L3")]).iloc[0]
        share = out["weights"].groupby("feature").size() / n_fits
        stable_here = set(share[share >= 0.8].index)
        rows.append({
            "drug": drug, "n_cell_lines": len(d.ids), "n_pten_altered": int(d.meta["pten_altered"].sum()),
            "r2_L1_pten": table.loc["L1", "r2_mean"], "r2_L3_pten_background": table.loc["L3", "r2_mean"],
            "r2_L4_expression": table.loc["L4", "r2_mean"],
            "gain_L4_over_L3": gain["gain_r2"], "gain_ci_low": gain["ci_low"], "gain_ci_high": gain["ci_high"],
            "stable_features": len(stable_here),
            "stable_features_shared_with_primary": len(stable_here & stable_primary),
            "share_of_primary_stable_features_also_stable": len(stable_here & stable_primary) / max(len(stable_primary), 1),
        })
    pd.DataFrame(rows).to_csv(args.results / "table_replication.csv", index=False)
    (args.results / "run_ladder_done.json").write_text(json.dumps({"drugs": drugs, "n_fits_per_drug": n_fits}))
    print("Ladder tables written", flush=True)


if __name__ == "__main__":
    main()
