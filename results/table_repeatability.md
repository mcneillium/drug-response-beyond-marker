EXPLORATORY: added after seeing the results of Experiments 1 to 6.
Agreement between two measurements of the same drug on the same cell lines (AUC).

| Drug | Repeat | Cell lines | Dose ranges (µM) | Pearson (95% interval) | Spearman | Best model R² here | Share of the ceiling |
| --- | --- | ---: | --- | --- | ---: | ---: | ---: |
| MK-2206 | GDSC1 against GDSC2 | 881 | 0.015625 to 4 and 0.010005 to 10 | 0.55 (0.48 to 0.62) | 0.44 | 0.150 | 27% |
| Uprosertib | within GDSC2 | 668 | 0.010005 to 10 and 0.010005 to 10 | 0.81 (0.76 to 0.85) | 0.73 | 0.180 | 22% |

"Best model R² here" is the elastic net with expression (L4) for that drug, from Experiment 5.
"Share of the ceiling" is that R² divided by the Pearson correlation.

| All repeated drugs | Number | Median Pearson | Lowest | Highest |
| --- | ---: | ---: | ---: | ---: |
| GDSC1 against GDSC2 | 122 | 0.45 | 0.01 | 0.79 |
| within GDSC2 | 9 | 0.78 | 0.07 | 0.95 |
| within GDSC2, same dose range only | 3 | 0.89 | 0.81 | 0.95 |
