128 PTEN-altered cell lines. k-means on 20 principal components of the 2,000 most variable genes.
Stability: mean adjusted Rand index (ARI) between the full clustering and 100 subsamples of 80%.
Reference: the same measure on data with no groups (mean and highest of 20 reference datasets).
Follows tissue: share of cell lines in the commonest group of their own tissue (tissues with 5 or more
cell lines), with the share expected if groups were unrelated to tissue in brackets.

| Version | k | Stability | Reference (mean, max) | Verdict | Group sizes | ARI with tissue | Follows tissue | ARI with centre | Response differs? p (AUC) | p (adjusted AUC) |
| --- | ---: | ---: | --- | --- | --- | ---: | --- | ---: | ---: | ---: |
| Expression as published | 2 | 0.66 | 0.87, 0.97 | in between | 108 / 20 | 0.08 | 1.00 (0.81) | 0.01 | 0.351 | 0.942 |
| Expression as published | 3 | 1.00 | 0.68, 0.84 | stable | 55 / 53 / 20 | 0.27 | 0.94 (0.49) | 0.00 | 3.08e-05 | 0.00921 |
| Expression as published | 4 | 0.93 | 0.62, 0.85 | stable | 51 / 36 / 21 / 20 | 0.37 | 0.91 (0.45) | 0.05 | 0.000745 | 0.0689 |
| Expression as published | 5 | 0.83 | 0.53, 0.66 | stable | 42 / 28 / 27 / 20 / 11 | 0.45 | 0.87 (0.37) | 0.01 | 9.5e-05 | 0.0262 |
| Tissue and centre removed | 2 | 0.75 | 0.80, 0.96 | in between | 103 / 25 | 0.06 | 0.91 (0.85) | 0.01 | 0.667 | 0.372 |
| Tissue and centre removed | 3 | 0.64 | 0.60, 0.85 | in between | 95 / 25 / 8 | 0.06 | 0.83 (0.77) | 0.01 | 0.18 | 0.261 |
| Tissue and centre removed | 4 | 0.40 | 0.51, 0.66 | unstable | 74 / 24 / 22 / 8 | 0.12 | 0.72 (0.60) | -0.00 | 0.000417 | 0.00286 |
| Tissue and centre removed | 5 | 0.57 | 0.44, 0.54 | in between | 59 / 34 / 15 / 12 / 8 | 0.07 | 0.55 (0.50) | 0.01 | 0.1 | 0.177 |
