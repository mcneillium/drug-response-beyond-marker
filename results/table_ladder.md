Ipatasertib, 924 cell lines. Held-out predictions, 5-fold cross-validation repeated 5 times.
Mean over repeats (standard deviation in brackets).

| Step | Inputs | R² | Spearman | Within-cancer-type Spearman |
| --- | --- | ---: | ---: | ---: |
| L0 | Training mean | -0.002 (0.000) | -0.042 (0.019) | -0.058 (0.064) |
| L1 | PTEN status | 0.051 (0.002) | 0.111 (0.014) | 0.074 (0.031) |
| L2 | Background factors | 0.068 (0.002) | 0.187 (0.012) | -0.249 (0.022) |
| L3 | PTEN + background | 0.112 (0.002) | 0.263 (0.012) | -0.073 (0.064) |
| L4 | + expression (elastic net) | 0.144 (0.013) | 0.282 (0.011) | 0.166 (0.032) |
| L5 | + expression (boosting) | 0.187 (0.010) | 0.321 (0.019) | 0.197 (0.067) |
| L1_strict | PTEN status, strict definition | 0.059 (0.003) | 0.100 (0.012) | 0.060 (0.042) |
| L3_strict | PTEN strict + background | 0.118 (0.004) | 0.265 (0.013) | -0.069 (0.056) |

**Gains in R² (95% bootstrap interval over cell lines)**

| Better | Simpler | Gain | 95% interval |
| --- | --- | ---: | --- |
| L1 | L0 | +0.053 | +0.019 to +0.090 |
| L3 | L2 | +0.043 | +0.016 to +0.071 |
| L3 | L1 | +0.061 | +0.017 to +0.099 |
| L4 | L3 | +0.033 | +0.005 to +0.060 |
| L5 | L3 | +0.076 | +0.014 to +0.130 |
| L5 | L4 | +0.043 | -0.003 to +0.086 |
| L4 | L1 | +0.093 | +0.053 to +0.129 |
| L5 | L1 | +0.137 | +0.059 to +0.204 |
