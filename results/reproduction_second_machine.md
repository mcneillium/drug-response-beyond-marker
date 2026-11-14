# Reproduction on a second machine, 14 November 2026

The tables in this folder come from the first run (Linux, 2 CPU cores,
Python 3.13), which a clean re-run on the same platform reproduced byte for
byte. On 14 November 2026 Paul McNeill re-ran every experiment on a second
machine (Windows, Python 3.14.7, NumPy 2.5.4, pandas 3.1.1, scikit-learn
1.9.1, LightGBM 4.7.0, SciPy 1.18.1) from commit `220d425`:
`run_ladder.py`, `run_overlap.py`, `run_clusters.py`,
`run_pathway_check.py`, `run_repeatability.py` and `plot_results.py`.

**Data.** All seven downloaded files have the same SHA-256 checksums as
`data_manifest.json`.

**Results.** Every result table has the same rows, columns and text. The
largest difference in any number:

| Table | Largest difference |
| --- | ---: |
| `table_clusters.csv`, `table_cluster_groups_by_tissue.csv`, `table_explanation_summary.csv`, `table_leak.csv` | 0 (identical) |
| `table_repeatability.csv` | 4e-16 |
| `table_feature_stability.csv`, `table_overlap.csv`, `table_ladder_gains.csv` | under 2e-9 |
| `table_pathway_check.csv`, `table_pathway_check_gains.csv` | under 4e-9 |
| `table_replication.csv` | 5e-7 |
| `table_ladder.csv` | 2e-6 |

The differences are floating-point rounding in the iterative elastic-net
fits. None changes a reported figure at the precision published, and no
conclusion changes. The first-run tables stay as the record.
