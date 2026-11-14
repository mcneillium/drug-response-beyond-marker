"""Experiment 3: stability of groups among PTEN-altered cell lines.

    python scripts/run_clusters.py

Two versions of the expression data are clustered: as published, and with the
average effect of tissue and sequencing centre removed (estimated from every
cell line in the study, not only the PTEN-altered ones).

Outputs (results/)
    table_clusters.csv / .md     stability, reference, drivers and response test per version and k
    raw/cluster_labels.csv       the full-data group of every cell line, per version and k
    raw/cluster_agreement.csv    the adjusted Rand index of every subsample
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal
from sklearn.metrics import adjusted_rand_score

from drbm import PRIMARY_DRUG, clusters, dataset

ROOT = Path(__file__).resolve().parents[1]
FACTORS = ["tissue_group", "source"]
VERSION_LABELS = {"as_published": "Expression as published", "adjusted": "Tissue and centre removed"}


def response_test(values: pd.Series, labels: np.ndarray) -> float:
    groups = [values[labels == g].to_numpy() for g in np.unique(labels)]
    groups = [g for g in groups if len(g) >= 3]
    return float(kruskal(*groups).pvalue) if len(groups) >= 2 else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "raw")
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    raw_dir = args.results / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    data = dataset.assemble(args.data, PRIMARY_DRUG)
    altered = data.meta.index[data.meta["pten_altered"]]
    meta = data.meta.loc[altered]

    # Genes that vary at all among the cell lines of the study (saves memory; changes nothing).
    varying = data.expression.columns[data.expression.var(axis=0) > 0]
    expression = data.expression[varying].astype("float64")
    versions = {
        "as_published": expression.loc[altered],
        "adjusted": clusters.remove_effects(expression, data.meta, FACTORS).loc[altered],
    }
    adjusted_auc = clusters.remove_effects(data.y.to_frame(), data.meta, FACTORS)["auc"].loc[altered]
    auc = data.y.loc[altered]

    rows, label_columns, agreements = [], {}, []
    for version, table in versions.items():
        full_points = clusters.components(table)
        for k in clusters.K_VALUES:
            labels, agreement = clusters.stability(
                lambda positions, table=table: clusters.components(table.iloc[positions]), len(table), k)
            reference = clusters.reference_stability(full_points, k)
            sizes = np.bincount(labels)
            rows.append({
                "version": version, "label": VERSION_LABELS[version], "k": k,
                "stability_mean_ari": agreement.mean(), "stability_sd": agreement.std(ddof=1),
                "stability_p05": np.percentile(agreement, 5),
                "reference_mean_ari": reference.mean(), "reference_max_ari": reference.max(),
                "clearly_above_reference": bool(agreement.mean() > reference.max()),
                "verdict": clusters.verdict(agreement.mean()),
                "smallest_group": int(sizes.min()), "group_sizes": " / ".join(str(s) for s in sorted(sizes, reverse=True)),
                "ari_with_tissue": adjusted_rand_score(meta["tissue"], labels),
                "share_in_commonest_group_of_own_tissue": clusters.tissue_agreement(meta["tissue"], labels)[0],
                "same_share_if_unrelated_to_tissue": clusters.tissue_agreement(meta["tissue"], labels)[1],
                "ari_with_centre": adjusted_rand_score(meta["source"], labels),
                "kruskal_p_auc": response_test(auc, labels),
                "kruskal_p_adjusted_auc": response_test(adjusted_auc, labels),
            })
            label_columns[f"{version}_k{k}"] = labels
            agreements.append(pd.DataFrame({"version": version, "k": k, "subsample": np.arange(len(agreement)),
                                            "ari": agreement}))
            print(f"  {version} k={k}: stability {agreement.mean():.2f}, reference {reference.mean():.2f}", flush=True)

    table = pd.DataFrame(rows)
    table.to_csv(args.results / "table_clusters.csv", index=False)
    pd.DataFrame(label_columns, index=altered).join(meta[["model_name", "tissue", "source"]]).to_csv(
        raw_dir / "cluster_labels.csv")
    pd.concat(agreements, ignore_index=True).to_csv(raw_dir / "cluster_agreement.csv", index=False)
    best = table[table["version"] == "as_published"].sort_values("stability_mean_ari").iloc[-1]
    pd.crosstab(meta["tissue"], pd.Series(label_columns[f"as_published_k{best['k']}"], index=altered, name="group")).to_csv(
        args.results / "table_cluster_groups_by_tissue.csv")

    lines = [f"{len(altered)} PTEN-altered cell lines. k-means on 20 principal components of the 2,000 most variable genes.",
             "Stability: mean adjusted Rand index (ARI) between the full clustering and 100 subsamples of 80%.",
             "Reference: the same measure on data with no groups (mean and highest of 20 reference datasets).",
             "Follows tissue: share of cell lines in the commonest group of their own tissue (tissues with 5 or more",
             "cell lines), with the share expected if groups were unrelated to tissue in brackets.", "",
             "| Version | k | Stability | Reference (mean, max) | Verdict | Group sizes | ARI with tissue | "
             "Follows tissue | ARI with centre | Response differs? p (AUC) | p (adjusted AUC) |",
             "| --- | ---: | ---: | --- | --- | --- | ---: | --- | ---: | ---: | ---: |"]
    for _, row in table.iterrows():
        lines.append(f"| {row['label']} | {row['k']} | {row['stability_mean_ari']:.2f} | "
                     f"{row['reference_mean_ari']:.2f}, {row['reference_max_ari']:.2f} | {row['verdict']} | "
                     f"{row['group_sizes']} | {row['ari_with_tissue']:.2f} | "
                     f"{row['share_in_commonest_group_of_own_tissue']:.2f} ({row['same_share_if_unrelated_to_tissue']:.2f}) | "
                     f"{row['ari_with_centre']:.2f} | "
                     f"{row['kruskal_p_auc']:.3g} | {row['kruskal_p_adjusted_auc']:.3g} |")
    (args.results / "table_clusters.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
