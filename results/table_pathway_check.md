EXPLORATORY: added after seeing the results of Experiments 1 to 5.
Ipatasertib, 924 cell lines, same folds. Mean over 5 repeats (standard deviation).
Pathway genes found in the expression data: 12 of 12.

| Step | Inputs | R² | Spearman | Within-cancer-type Spearman |
| --- | --- | ---: | ---: | ---: |
| L1 | PTEN status (from Experiment 1) | 0.051 (0.002) | 0.111 (0.014) | 0.074 (0.031) |
| P1 | PTEN gene expression level alone | 0.000 (0.003) | 0.037 (0.006) | 0.120 (0.037) |
| P2 | PTEN status + PTEN expression level | 0.045 (0.004) | 0.100 (0.013) | 0.113 (0.038) |
| L3 | PTEN status + background (from Experiment 1) | 0.112 (0.002) | 0.263 (0.012) | -0.073 (0.064) |
| P3 | PTEN status + background + 12 pathway genes | 0.132 (0.004) | 0.255 (0.009) | 0.116 (0.029) |

| Better | Simpler | Gain in R² | 95% interval |
| --- | --- | ---: | --- |
| P1 | L1 | -0.050 | -0.084 to -0.017 |
| P2 | L1 | -0.006 | -0.010 to -0.001 |
| P3 | L3 | +0.021 | +0.004 to +0.037 |
