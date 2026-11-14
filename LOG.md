# Log

Every decision and every run, including what went wrong. Newest section last.

---

## 2026-11-13: how this study was run

Paul asked for the portfolio studies to be built without stopping for
approval. Everything below was done in one unattended
session, in a cloud Linux container (2 CPU cores, 7 GB memory, no GPU,
Python 3.13.18). Nothing here has yet been re-run on Paul's own machine
(a clean re-run in the same kind of environment is recorded below).

**Every choice in the tables below is for Paul to confirm or change.**

## Data check (before the protocol was written)

- The plan named the Broad Institute's DepMap portal. Its download page sits
  behind an interactive browser check and could not be scripted. That check
  was not worked around. The Sanger sources were used instead: GDSC2 release
  8.5 for drug response, Cell Model Passports for everything else.
- Four AKT inhibitors are in the screen. After the filters below, 867 to 934
  cell lines are usable per drug.
- Only 6 prostate cell lines were screened with Ipatasertib. A prostate-only
  study is not possible with this data, so all cancer types are used.
- Gene expression came from two sequencing centres (Broad 535 cell lines,
  Sanger 389 among those used). Treated as a technical factor throughout.
- Protein data covers 901 of the 924 cell lines used. Partial overlap
  therefore had to be simulated by hiding protein data.
- The expression file holds a small number of cells with two numbers in them
  (for example `0.0144 0.0704`). They are treated as missing, and any gene
  with a missing value is dropped. 36,423 genes remain.
- AUC is strongly skewed: median 0.93, and two thirds of cell lines are above
  0.90 (barely affected by the drug). PTEN-altered lines have a median of
  0.89 against 0.94 for the rest.

## Choices made during the run

| # | Choice | Reason | Alternative |
| --- | --- | --- | --- |
| D1 | Sanger Cell Model Passports and GDSC in place of the DepMap portal | The portal could not be downloaded by script | Download DepMap files by hand |
| D2 | Ipatasertib as the primary drug | An AKT inhibitor trialled in prostate cancer with PTEN loss (IPATential150, Sweeney et al., Lancet 2021; checked 13 Nov against the journal record) | MK-2206, the most studied AKT inhibitor in cell lines |
| D3 | MK-2206, Uprosertib and Afuresertib for replication | The other AKT inhibitors in the same screen | A truly independent dataset (CTRP or PRISM); not done |
| D4 | AUC as the response, not IC50 | IC50 is extrapolated for cell lines the drug barely affects | ln(IC50) |
| D5 | PTEN "altered" = a cancer-driver mutation or both copies deleted | Uses both kinds of evidence in the data | The stricter definition (D6) as primary |
| D6 | Strict definition (truncating mutation or deletion) as a sensitivity check only | Declared in advance | See result note 2 |
| D7 | A cell line with no mutation or copy number data has no PTEN status and is left out | Missing is not the same as normal | Assume normal |
| D8 | Where a cell line was sequenced by both centres, the Sanger run is used | Needed a rule; applied before any modelling | Broad run, or average |
| D9 | Background factors = cancer type (types under 10 cell lines grouped as "Other") + sequencing centre | The two obvious confounders | Add growth medium, microsatellite status |
| D10 | Folds: 5 x 5, stratified by tissue (tissues under 5 cell lines grouped), seeds 0 to 4 | Standard; the same cell lines are held out for every model | Leave-one-tissue-out, which asks a harder question |
| D11 | The 2,000 most variable genes, chosen inside each training fold | A filter that does not look at the response | More genes, or a knowledge-based list |
| D12 | Elastic net mixing 0.5; strength by inner 5-fold cross-validation | Middle of the range; declared in advance | Tune the mixing too |
| D13 | Elastic net speed settings (30 strengths over a factor of 100, tolerance 0.001, random coordinate order) | 77 s per fit became 3 s, with the same held-out R² on the fold used for timing (0.170 and 0.171) | Library defaults, about 25 times slower |
| D14 | Gradient boosting settings fixed and not tuned (400 trees, learning rate 0.03, 7 leaves) | A baseline should be plain | Tune inside the folds |
| D15 | Leak demonstration: 200 genes by correlation, ridge regression, one shuffle of the response | Simple and common in practice | Several shuffles; see limits |
| D16 | Partial overlap simulated at 97%, 75%, 50% and 25% by hiding protein at random | The real overlap is 97% | Hide by tissue or by batch, which is more realistic |
| D17 | Gaps inside the protein data (after the 80% filter) filled with the training median | Needed for any linear model | A proper imputation model |
| D18 | Late-fusion weight from 0, 0.1, ..., 1 by inner 5-fold cross-validation | Keeps the choice inside the training fold | A stacked model |
| D19 | Clustering: k-means on 20 principal components of 2,000 genes, k = 2 to 5 | Plain and common | Hierarchical or consensus clustering |
| D20 | "Recluster" repeats every step on each subsample | Stricter than reclustering fixed components | Fixed components |
| D21 | Stability thresholds 0.75 (stable) and 0.50 (unstable) | A judgement, fixed in advance | Any other |
| D22 | A no-groups reference for stability (amendment 6) | k-means gives repeatable splits of structureless data | None; see result note 6 |
| D23 | Tissue and centre effects removed by least squares fitted on all 924 cell lines | Uses the most data to estimate each effect | Fit on the 128 only |
| D24 | Top-20 features by mean absolute SHAP on the held-out fold | Declared in advance | Top 10 or top 50 |
| D25 | Pathway gene list from general knowledge | A sanity check | A curated pathway database |
| D26 | Gains computed within each repeat, then averaged (amendment 9) | Corrects an error; see "what went wrong" | None |
| D27 | Noise ceiling read as the Pearson correlation between two measurements (amendment 10) | Standard reading when both have equal, independent noise | Report only the raw correlation and draw no ceiling |

## Runs

| # | Script | What it did | Time |
| --- | --- | --- | --- |
| 1 | `run_ladder.py` | Ladder for 4 drugs (100 elastic-net fits, 25 boosting fits), leak, stability, replication | about 8 min alone (20 min in the first run, alongside other work) |
| 2 | `run_overlap.py` | 4 shares x 25 folds x 4 strategies | about 6 min |
| 3 | `run_clusters.py` | 2 versions x 4 values of k x 100 subsamples, plus 20 reference datasets each | about 2 min |
| 4 | `run_pathway_check.py` | Exploratory: pathway genes forced in | under 1 min |
| 5 | `run_repeatability.py` | Exploratory: agreement between repeat measurements (131 drug pairs) | about 2 min, mostly reading files |
| 6 | `plot_results.py` | Five figures | seconds |

Tests: 54, all passing. They use small hand-made data and need no download.

## Results (details in docs/report.md)

1. **Ladder, Ipatasertib.** R²: PTEN 0.051, background 0.068, PTEN +
   background 0.112, + expression with elastic net 0.144, + expression with
   boosting 0.187. Gain from expression over PTEN + background: +0.033 (95%
   interval +0.005 to +0.060) with elastic net, +0.076 (+0.014 to +0.130) with
   boosting. Boosting over elastic net: +0.043 (-0.003 to +0.086), not clear.
2. **Strict PTEN definition.** Slightly better than the primary one (R² 0.059
   against 0.051 alone; 0.118 against 0.112 with background). It was declared
   as a sensitivity check, so the primary definition was not switched.
3. **Leak.** Genes chosen with all cell lines: R² 0.210 on the real response
   and 0.137 on a shuffled response where nothing is real. Chosen inside the
   folds: -0.005 and -0.173.
4. **Partial overlap.** Protein data added nothing. Late fusion matched
   expression only (0.152 to 0.161 against 0.156). Complete cases and
   fill-the-gaps were worse at every share, and collapsed at 25% (R² about
   0.03 on the cell lines with protein).
5. **Clusters as published.** k = 3 is perfectly stable (ARI 1.00; reference
   0.68, highest 0.84). The groups are blood (20), a group led by brain and
   skin (53), and a group of carcinomas (55). 94% of cell lines sit in the
   commonest group of their own tissue.
6. **Clusters after removing tissue and centre.** k = 2: 0.75 against a
   reference of 0.80. k = 3: 0.64 against 0.60. k = 4: 0.40 against 0.51.
   k = 5: 0.57 against 0.44 (highest reference 0.54). Nothing is both stable
   and clearly above the reference. Note that the reference itself reaches
   0.80 to 0.87 at k = 2: data with no groups looks "stable" there.
7. **Explanations.** Top-20 lists from 25 boosting fits overlap by 0.20 on
   average. Two features are in the top 20 in at least 80% of fits: PTEN
   status (all 25) and the gene EMILIN2 (22). 134 different features reach a
   top-20 list at least once. The elastic net selects 54 features per fit on
   average; 19 are selected in at least 80% of fits.
8. **Pathway sanity check.** 10 of the 12 declared pathway genes, including
   PTEN itself, were never model inputs: they are not among the 2,000 most
   variable genes. The check was close to empty as designed.
9. **Replication.** Gain from expression (elastic net over PTEN +
   background): Ipatasertib +0.033 (+0.005 to +0.060), MK-2206 +0.008
   (-0.021 to +0.038), Uprosertib +0.082 (+0.045 to +0.118), Afuresertib
   +0.015 (-0.011 to +0.043). Of Ipatasertib's 19 most often selected
   features, 5, 12 and 8 are also among the most often selected for the
   others.
10. **Exploratory pathway check.** PTEN's own expression level alone: R²
    0.000. PTEN status + background + 12 pathway genes: 0.132 against 0.112
    (gain +0.021, +0.004 to +0.037).

11. **Exploratory repeatability.** Uprosertib, measured twice within GDSC2 at
    the same doses on 668 shared cell lines: Pearson 0.81 (0.76 to 0.85). The
    best model for that drug explains 0.180, about 22% of that ceiling.
    MK-2206, GDSC1 against GDSC2 (different assay and doses), 881 cell lines:
    0.55 (0.48 to 0.62). Across all repeated drugs: median 0.78 within GDSC2
    (9 drugs, range 0.07 to 0.95; 0.89 for the 3 with the same dose range),
    median 0.45 between screens (122 drugs, range 0.01 to 0.79).

## Against the negative-result criteria fixed in the protocol

| Criterion | Outcome |
| --- | --- |
| Expression gain interval includes zero | Not met for Ipatasertib (both models exclude zero, the elastic net only just). Met for 2 of the 3 replication drugs |
| No clustering reaches ARI 0.50 after adjustment | Not met as written: three values of k exceed 0.50. But the criterion had no reference. Against the no-groups reference added before the run, none is both stable and clearly above it. I read this as "no trustworthy subgroups", and say so with this caveat |
| Top-20 overlap below 0.3 | Met: 0.20. The explanations are not stable enough to hand to a biologist |

## What went wrong, and what was changed

- **Gains were first computed the wrong way.** Predictions were averaged over
  the 5 repeats before comparing R². That is a small ensemble and flattered
  the less stable models (boosting over L3: +0.094 instead of +0.076). Found
  when the gain did not match the table. Fixed (amendment 9), with a test that
  shows the size of the old error. All gains in the results use the fixed method.
- **Feature importance was first averaged only over the fits where a gene was
  an input.** A gene that entered 3 fits out of 25 ranked second overall.
  Fixed: importance is averaged over all 25 fits, counting zero where the gene
  was not an input. The fix did not change the top-20 overlap (0.20), which is
  computed per fit.
- **The overlap experiment was too slow** (about 2 minutes per fold) because
  of column-by-column table operations. Rewritten with plain arrays; the
  expression half of late fusion is now computed once per fold. No results
  had been read from the slow version; it was stopped after 3 folds.
- **The adjusted Rand index against tissue was a poor measure** (amendment
  7). It gave 0.27 for groups that follow tissue almost exactly.
- **The pathway sanity check was nearly empty** (result 8). This is a design
  flaw in the protocol, reported as such. Experiment 6 was added to see what
  it cost, and is labelled exploratory.

## Reproducibility checks

- Two folds of the ladder (repeat 0 fold 0, repeat 3 fold 2) were refitted
  from scratch: elastic-net and boosting predictions matched the saved ones
  exactly.
- The leak experiment and the cluster experiment were each run twice in full
  and gave identical tables.
- The overlap experiment was run once in its final form. Not repeated.
- The full ladder for all four drugs was first run once.
- **Clean re-run, 13 November 2026 (evening).** The committed repository was
  cloned into an empty folder. A new Python environment was built from
  `requirements-lock.txt`, the data was downloaded again, and every script
  was run in the order given in the README (about 18 minutes in total).
  - The checksums of all seven downloaded files matched `data_manifest.json`.
  - All 54 tests passed.
  - Every table in `results/` came out byte-for-byte identical.
  - The five figures show the same numbers in a different typeface: the
    first environment had the Inter font configured and a clean one does not.
  - This was the same machine type, operating system and Python version as
    the first run. It shows the repository is complete and the results are
    deterministic. It is not a re-run on a second machine.

## Limits

- **Cell lines, not patients.** Nothing here says anything about treatment.
- **One screen.** The replication drugs come from the same screen and the
  same cell lines. Their responses correlate with Ipatasertib's at 0.75 to
  0.89. Agreement shows consistency, not independent confirmation.
- **The intervals understate uncertainty.** They resample cell lines. They do
  not include the change from refitting a model, and there are only 924 cell
  lines.
- **R² is driven by a minority.** Two thirds of cell lines barely respond.
  A few badly missed sensitive cell lines can move R² a lot (this is what
  happened to the correctly done ridge model in the leak experiment).
- **Partial overlap was simulated at random.** In real studies the missing
  samples differ systematically from the rest.
- **The noise ceiling is by analogy.** Ipatasertib was measured once. The
  ceiling comes from a related drug, and rests on the two measurements having
  equal and independent noise. Repeats within one laboratory can share errors,
  which would make the ceiling too high.
- **One shuffle** in the leak demonstration.
- **No multiple-testing correction** on the 16 Kruskal-Wallis tests in the
  cluster table. None of them is used for a claim.
- **Boosting and elastic net were not tuned.** A tuned model might explain
  more. It would not change the instability findings by itself.
- **No biological interpretation** has been checked by a biologist. Gene
  names are reported as what the models used, nothing more.

## Not verified

- The terms of use for the GDSC and Cell Model Passports data (CITATION.md).
- The citations for the merged expression data.
- The papers in CITATION.md were not read in full; only their journal details
  were checked.
- The pathway gene list.

## 2026-11-13 (evening): question, authorship statement and first report added

At Paul's instruction ("you do everything") the research question, a
status note and
`docs/report.md`. The README says no prediction was recorded in advance.
The report's related-work section is left unwritten because the papers have
not been read in full.

## 2026-11-14: prepared for publication

- The names of individual academics and references to a specific project
  advert were removed from `docs/protocol.md`. No design content changed.
- Personal study notes (the walkthrough) and the working rules file are
  kept out of the public repository, under `notes/`, which git ignores.
- The public history starts from a single commit. The earlier commits are
  kept privately.
- Data licence terms were checked and recorded in `CITATION.md`.

## 2026-11-14: corrections from a review of the logged choices

- The report said "the only stable subgroups are tissue". That claims more
  than one clustering method can show. It now says no stable subgroups beyond
  tissue were found with these methods.
- Gene selection by variance and the standardisation are fitted on each outer
  training fold, and the elastic net then chooses its strength by inner
  cross-validation on that matrix. The inner validation rows therefore helped
  choose the 2,000 genes. This cannot leak the outer test fold, so the
  held-out R² values stand, but the tuning is not fully nested. Recorded under
  D11 for Paul's decision.
- The leak demonstration uses one shuffle of the response. Repeated
  permutations would give a null distribution. Recorded under D15 for Paul's
  decision.

## 2026-11-14: reproduced on a second machine

Paul re-ran all six scripts on his Windows machine (Python 3.14.7) from
commit `220d425`. The downloaded data matched the manifest checksums. Every
result table matched the published one to within 2e-6, and four were
identical. Details: `results/reproduction_second_machine.md`.

