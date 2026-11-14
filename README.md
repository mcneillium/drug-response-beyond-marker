# Drug response beyond a single marker

**Status:** all seven experiments have been run end to end on public data,
and a clean re-run from a fresh clone reproduced every result table exactly
(same kind of machine; see `LOG.md`). A first report is in `docs/report.md`.
Every table was reproduced on a second machine (Windows, Python 3.14) on
14 November 2026: see `results/reproduction_second_machine.md`.

## Question

For a drug where the status of one gene (PTEN) is a known but incomplete
guide to response, how much more of the response does wider molecular data
explain, and are the patterns found stable enough to trust?

<!-- TODO (Paul): rewrite the question in your own words. -->

## What I expected to find

No prediction was written down before the results existed, so none is claimed
here. The protocol did fix, in advance, what would count as a negative result;
`LOG.md` reports the outcome against each of those criteria.

## The setup in brief

- **Setting.** One gene, PTEN, is a known but incomplete guide to how cancer
  cells respond to drugs that block the AKT protein. This study measures how
  much more of the response wider molecular data explains, and whether the
  patterns found can be trusted.
- **Data.** 924 cancer cell lines from the GDSC drug screen with molecular
  data from the Sanger Cell Model Passports. All public. See `CITATION.md`.
- **Drug.** Ipatasertib (primary). Three other AKT inhibitors are used to see
  whether findings repeat: MK-2206, Uprosertib, Afuresertib.
- **Response.** AUC, the area under the dose-response curve. 1 means no
  effect; lower means more sensitive.
- **Marker.** PTEN "altered" (a cancer-driver mutation or both copies
  deleted): 128 of the 924 cell lines.
- **Honesty rules.**
  - The protocol (`docs/protocol.md`) was written before any model was fitted.
    Later changes are listed in it as dated amendments, and those made after
    seeing results are labelled exploratory.
  - Every step that learns from data (choosing genes, scaling, filling gaps,
    tuning) happens inside the training folds only.
  - The same saved folds are used for every model: 5-fold cross-validation
    repeated 5 times, in `splits/folds_v1.json`.

## The experiments

| # | Question | Script |
| --- | --- | --- |
| 1 | How much does each kind of information add? PTEN, then cancer type and sequencing centre, then gene expression | `run_ladder.py` |
| 1b | How large is the error from choosing genes before cross-validating? | `run_ladder.py` |
| 2 | Protein data exists for only some cell lines. Which way of handling that works? (simulated) | `run_overlap.py` |
| 3 | Are there stable subgroups among PTEN-altered cell lines, or only tissue? | `run_clusters.py` |
| 4 | Do the models point to the same genes each time they are refitted? | `run_ladder.py` |
| 5 | Does the main result repeat on three other AKT inhibitors? | `run_ladder.py` |
| 6 | Exploratory, added after seeing results: do the known pathway genes help if forced in? | `run_pathway_check.py` |
| 7 | Exploratory, added last: how repeatable is the drug response itself, where it was measured twice? | `run_repeatability.py` |

## Results in one screen

All numbers are held-out. R² is the share of the variation in response that a
model explains (1 = all of it, 0 = none).

| Finding | Number |
| --- | --- |
| PTEN status alone | R² 0.05 |
| PTEN + cancer type + sequencing centre | R² 0.11 |
| Adding 2,000 genes, elastic net | R² 0.14 (gain +0.033, 95% interval +0.005 to +0.060) |
| Adding 2,000 genes, gradient boosting | R² 0.19 (gain +0.076, 95% interval +0.014 to +0.130) |
| Choosing genes before cross-validating, on shuffled data with nothing to find | R² 0.14, entirely false |
| Protein data on top of expression | no gain at any level of overlap |
| Stable groups among PTEN-altered lines | only before tissue is removed, and 94% of cell lines sit with their own tissue |
| Overlap between the top-20 gene lists of refitted models | 0.20 (my pre-set bar for "stable enough" was 0.30) |
| Expression gain on the three other AKT inhibitors | interval excludes zero for 1 of 3 |
| Agreement between two measurements of a related drug (Uprosertib) | 0.81, so about 80% is repeatable and the model explains 18% |

So the short version is mostly negative: wider data helps a little, more than
four fifths of the response stays unexplained (and noise in the measurements
accounts for only a small part of that), the subgroups are tissue, and the
gene lists are not stable. `docs/report.md` sets out each result and
its limits. Figures are in `figures/`.

## Reproduce

```
conda env create -f environment.yml
conda activate drbm
pip install -e .[dev,analysis]

python -m drbm.download --dest data/raw   # about 1.4 GB, no login needed
pytest                                    # 54 tests, no data needed

python scripts/run_ladder.py              # about 8 minutes on 2 cores
python scripts/run_overlap.py             # about 6 minutes
python scripts/run_clusters.py            # about 2 minutes
python scripts/run_pathway_check.py       # under a minute
python scripts/run_repeatability.py       # about 2 minutes
python scripts/plot_results.py
```

The first run reads a 1.2 GB archive and needs about 6 GB of memory for a few
minutes; later runs use a cached copy. Run the scripts one at a time on a
machine with 8 GB. `data_manifest.json` holds the checksums of the files the
results here were produced from; compare it with your `data/raw/MANIFEST.json`.

## Where things are

```
docs/protocol.md         what was fixed in advance, and every later amendment
docs/report.md           the study in one document
LOG.md                   every choice and every run, including what went wrong
splits/folds_v1.json     the fixed folds for all four drugs, with a fingerprint
src/drbm/download.py     fetch the data and record checksums
src/drbm/load.py         read the raw files into tidy tables
src/drbm/dataset.py      assemble one drug's tables
src/drbm/folds.py        make, save and load the folds
src/drbm/scores.py       R2, rank correlations, bootstrap intervals
src/drbm/models.py       the baseline ladder and the leak demonstration
src/drbm/overlap.py      four ways to handle a partly missing data type
src/drbm/clusters.py     cluster stability, the no-groups reference, tissue checks
src/drbm/repeats.py      agreement between repeat measurements (noise ceiling)
scripts/                 one script per experiment, plus the figures
results/                 the tables (large intermediate files are not committed)
figures/                 five figures
tests/                   54 tests on small hand-made data
```

## Status of this work

Author: Paul McNeill. This is work in progress and has not been peer reviewed.
`LOG.md` lists the 27 choices made during the unattended runs (D1 to D27),
for me to confirm or change, and the mistakes found and fixed. As of
13 November 2026 I have not yet re-run the study on my own machine or
confirmed those choices.

<!-- TODO (Paul): restore a statement of how any AI tools were used before this goes to anyone else, and keep it accurate. -->

## Licence

Code: MIT (see `LICENSE`). Data: not included; see `CITATION.md`.
