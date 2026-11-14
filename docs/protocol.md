# Protocol: drug response beyond a single marker

Written on 13 November 2026 **before** any model was fitted. Only the data
checks in section 1 had been run. Its purpose is to fix the choices in
advance, so that nothing is tuned on the results.

## The question

For a drug where the status of one gene (PTEN) is a known but incomplete
guide to response, how much more of the response does wider molecular data
explain, and are the patterns found stable enough to trust?

This mirrors, on public data, a common setting in treatment-resistance
research: response that one marker does not fully explain, several types of measurement
that do not cover every sample, and a warning that groups found by clustering
may be unstable or driven by technical differences.

## 1. Data (checked before this was written)

| Item | Source | What the check found |
| --- | --- | --- |
| Drug response | GDSC2, release 8.5 | 286 drugs, 969 cell lines |
| Cell line details, mutations, copy number, gene expression, protein abundance | Sanger Cell Model Passports | Same cell line identifiers as GDSC |

- **Drug.** Ipatasertib, an AKT inhibitor, screened on 958 cell lines. AKT sits
  directly downstream of PTEN, so PTEN loss is a recognised marker for this
  drug class. Ipatasertib was chosen over the other AKT inhibitors in the
  screen because it has been trialled in prostate cancer with PTEN loss (my
  recollection; to be verified by Paul). The choice was made before any
  modelling.
- **Replication drugs.** MK-2206, Uprosertib and Afuresertib, three other AKT
  inhibitors in the same screen. Their responses correlate with Ipatasertib's
  at 0.75 to 0.89, so agreement across them shows consistency within one
  screen. It is not an independent dataset.
- **Response measure.** AUC, the area under the dose-response curve: 1 means
  the drug had no effect, lower means more sensitive. It is used in preference
  to IC50, which is extrapolated for cell lines the drug barely affects.
- **Cell lines used.** The 924 with a response, a known PTEN status and gene
  expression. All cancer types are used: only 6 prostate cell lines were
  screened, far too few to model.
- **PTEN status.** "Altered" means any cancer-driver mutation in PTEN or both
  copies deleted: 128 of 924 cell lines. A stricter definition (truncating
  mutation or deletion, 99 cell lines) is a sensitivity check. A cell line with
  no mutation or copy number data has no status; it is not assumed normal.
- **Two technical facts to control for.**
  - Cancer type predicts drug response by itself. Rare cancer types (fewer
    than 10 cell lines) are grouped as "Other".
  - Gene expression came from two sequencing centres (Broad for 535 cell
    lines, Sanger for 389). This is a technical difference of exactly the kind
    that can pass for biology.
  These two together are called the **background factors**.
- **Protein data.** 8,457 proteins, 39% of values missing. It covers 901 of
  the 924 cell lines, so the two data types overlap almost completely. Partial
  overlap is therefore studied by hiding the protein data for a random share
  of cell lines (section 4), and this is always labelled as simulated.

## 2. How everything is scored

- **Folds.** 5-fold cross-validation repeated 5 times (seeds 0 to 4),
  stratified by tissue so every fold has a similar mix. The folds are saved
  once in `splits/folds_v1.json`.
- **Every data-dependent step happens inside the training folds**: choosing
  which genes to use, scaling, filling gaps, tuning. The held-out fold is used
  only to predict.
- **Scores**, all computed on held-out predictions:
  - R²: the share of the variation in response explained. This is the main
    score, because the question asks "how much".
  - Spearman rank correlation.
  - Within-cancer-type Spearman: the correlation inside each cancer type with
    at least 20 cell lines, then the median over types. This shows whether a
    model does more than tell cancer types apart.
- **Uncertainty.** Mean and standard deviation of each score over the 5
  repeats, and a 95% bootstrap interval (2,000 resamples of cell lines) for
  each gain over a simpler model.

## 3. Experiment 1: the baseline ladder

Each step adds one thing. All use the same folds.

| Step | Inputs | Model |
| --- | --- | --- |
| L0 | none | the training-fold mean |
| L1 | PTEN status | two group means |
| L2 | background factors | ridge regression |
| L3 | PTEN status + background | ridge regression |
| L4 | L3 + gene expression | elastic net |
| L5 | L3 + gene expression | gradient boosting (LightGBM) |

- **Genes.** The 2,000 genes whose expression varies most across the training
  cell lines. This filter does not look at the drug response.
- **Elastic net.** Mixing parameter 0.5; strength chosen by 5-fold
  cross-validation inside the training fold.
- **Gradient boosting.** Fixed, plain settings, not tuned (in `models.py`).
- **The headline number** is the gain in R² from L3 to L4 or L5: what
  expression adds beyond PTEN and the background factors.

### Experiment 1b: the leak, measured

The most common error with thousands of features and hundreds of samples is to
choose features using all the data before cross-validating. To show its size:

- Pick the 200 genes most correlated with response, (a) wrongly, using every
  cell line once at the start, and (b) correctly, inside each training fold.
- Fit ridge regression on those genes and score both ways.
- Repeat with the response values shuffled between cell lines, so there is
  nothing real to find. The correct method should score about zero. Whatever
  the wrong method scores on shuffled data is pure illusion.

## 4. Experiment 2: when a data type is missing for some cell lines

Protein data is kept for a share of cell lines and hidden for the rest:
0.97 (as published), 0.75, 0.50 and 0.25, chosen at random with the repeat's
seed. Four strategies are compared, all with ridge regression (strength chosen
inside the training fold) so that only the strategy differs:

| Strategy | Trains on | Can predict |
| --- | --- | --- |
| Expression only | all training cell lines | every cell line |
| Complete cases | training cell lines with both data types | only cell lines with both |
| Fill the gaps | all training cell lines; missing protein values set to the training median, with a flag | every cell line |
| Late fusion | one model per data type, each on the cell lines that have it; predictions blended with a weight chosen inside the training fold | every cell line; expression model alone where protein is missing |

Proteins are used if detected in at least 80% of the training cell lines that
have protein data. Reported for each share: R² on all cell lines, on those
with both data types, and on those without protein, and how many cell lines
each strategy could use.

## 5. Experiment 3: are there stable groups among PTEN-altered cell lines?

- **Data.** The 128 PTEN-altered cell lines; the 2,000 most variable genes
  among them; scaled; reduced to 20 principal components.
- **Two versions.** (a) Expression as it is. (b) Expression with the average
  effect of tissue and sequencing centre removed, estimated from all 924 cell
  lines. Version (a) is expected to recover tissue; version (b) asks whether
  anything is left.
- **Clustering.** k-means with k = 2, 3, 4, 5.
- **Stability.** Recluster 100 random subsamples of 80% of the cell lines and
  compare each with the full clustering on the cell lines they share, by
  adjusted Rand index (ARI: 1 = identical groups, 0 = no better than chance).
  My working convention, fixed now: a mean ARI above 0.75 counts as stable,
  below 0.50 as unstable. These thresholds are a judgement, not a standard.
- **What drives the groups.** ARI between the groups and tissue, and between
  the groups and sequencing centre.
- **Do the groups differ in response?** Kruskal-Wallis test on AUC across
  groups, for the most stable k in each version, and again on AUC after
  removing the tissue and centre effects.

## 6. Experiment 4: are the explanations stable?

- From the 25 elastic-net fits of L4: how often each gene is given a non-zero
  weight.
- From the 25 boosting fits of L5: each feature's mean absolute SHAP value on
  the held-out fold; the top 20 per fit; how often each feature is in the top
  20; and the average overlap (Jaccard index) between the top-20 lists of
  different fits.
- **A sanity check, declared now.** Where do PTEN status and these pathway
  genes rank: PTEN, AKT1, AKT2, AKT3, PIK3CA, PIK3CB, PIK3R1, MTOR, TSC1, TSC2,
  INPP4B, PDPK1? The list is from general knowledge of the PI3K/AKT pathway and
  should be checked by a biologist.
- No biological claim will be made. "The model relies on these genes" is the
  strongest statement allowed.

## 7. Experiment 5: does it repeat on the other AKT inhibitors?

The ladder steps L1, L3 and L4 are repeated for the three replication drugs.
Reported: the gain from expression for each, and how many of Ipatasertib's
most often selected genes are also selected for the others.

## What would count as a negative result

- The gain from L3 to L4 and L5 has a 95% interval that includes zero: then,
  by these methods, expression adds nothing beyond PTEN and the background.
- No clustering reaches a mean ARI of 0.50 after removing tissue and centre:
  then there are no trustworthy subgroups here.
- Top-20 lists overlap by less than 0.3 on average: then the explanations are
  not stable enough to hand to a biologist.

Any of these would be reported as it is.

## Amendments

Each is dated and says what had and had not been run when it was made.

### 13 November 2026, before Experiments 2 and 3 were run

Experiment 1 was running for the first time; no result table had been read.

1. **Elastic net settings (speed only).** One fit took 77 seconds with the
   library defaults, which would have made the ladder take hours. The fit now
   uses 30 strengths spanning a factor of 100, a looser stopping tolerance and
   random coordinate order. On the one fold used for timing, held-out R² was
   unchanged (0.170 and 0.171). The mixing parameter (0.5) and the inner
   5-fold selection are as written above.
2. **Experiment 2, gaps inside the protein data.** Even a cell line that has
   protein data lacks some proteins. After the 80% filter those remaining gaps
   are filled with the training-fold median of that protein, in every strategy
   that uses protein. The "fill the gaps" strategy differs only in also
   filling whole cell lines that have no protein data, and adding a flag.
3. **Experiment 2, ridge strengths.** With several thousand inputs the best
   strength can be larger than in Experiment 1, so the grid is 1 to 1,000,000
   (13 values). The strength chosen in each fit is saved, to check it is not
   stuck at either end.
4. **Experiment 2, the blend weight.** Chosen from 0, 0.1, ..., 1 by 5-fold
   cross-validation inside the training fold, using the training cell lines
   that have protein data.
5. **Experiment 3, what "recluster" means.** Each subsample repeats every
   step: choosing the 2,000 genes, scaling, principal components, k-means.
   This is stricter than reclustering fixed components.
6. **Experiment 3, a reference for "stable".** k-means can give repeatable
   groups on data that has no groups at all, simply by cutting along the
   direction of greatest spread. So the same stability measure is also
   computed on reference data with no groups: random draws from one bell-shaped
   cloud with the same spread along each of the 20 components (20 such
   reference datasets per k). A clustering only counts as evidence of real
   groups if it is clearly more stable than this reference. The smallest group
   size is also reported, because a "group" of one or two cell lines can be
   perfectly stable and still mean nothing.

### 13 November 2026, after Experiments 1 to 5 had been run and read

These two additions were made **after seeing results**. They are labelled
exploratory wherever they are reported and no headline claim rests on them.

7. **Experiment 3, a second measure of "follows tissue".** The adjusted Rand
   index against tissue turned out to be a poor measure: three groups that
   each contain several whole tissues score low on it although they follow
   tissue almost exactly. Added: the share of cell lines that sit in the
   commonest group of their own tissue (tissues with at least 5 cell lines),
   beside the share expected if groups were unrelated to tissue (labels
   shuffled 200 times). The original measure is still reported.
8. **Experiment 6 (exploratory): the pathway genes, forced in.** The sanity
   check in Experiment 4 was largely empty: 10 of the 12 declared pathway
   genes, including PTEN itself, are not among the 2,000 most variable genes,
   so the models never saw them. To see what that filter cost, on the same
   folds: (a) the expression level of the PTEN gene alone, in place of PTEN
   mutation status; (b) PTEN status + background + the 12 pathway genes, with
   ridge regression. Compared with L1 and L3.
9. **A correction to how gains are computed.** The first version of the
   results averaged each model's held-out predictions over the 5 repeats and
   then compared R². An average of five predictions is a small ensemble, so
   this flattered the less stable models (it gave boosting a gain of +0.094
   over L3 where the difference between the two mean R² values in the table
   is +0.075). Gains are now computed within each repeat and then averaged,
   so they agree with the table. The interval still resamples cell lines
   (2,000 times), using the same resample for every repeat. The first version
   is described here so the change is visible; its numbers are not used.

### 13 November 2026, later: one more exploratory experiment

10. **Experiment 7 (exploratory): how repeatable is the response itself?**
    Added after all results above were known, prompted by the finding that
    more than 80% of the response is unexplained. Some of that must be
    measurement noise, which no model can explain. Where the same drug was
    measured twice on the same cell lines, the agreement between the two
    measurements estimates how much of the variation is real. Two kinds of
    repeat exist in the public data:
    - within the GDSC2 screen: 9 drugs were screened under two drug ids,
      including the AKT inhibitor Uprosertib;
    - across screens: 122 drugs are in both GDSC1 (an older screen with a
      different assay) and GDSC2, including the AKT inhibitor MK-2206.
    For each pair: the cell lines measured both times, and the Pearson and
    Spearman correlation of AUC. Reading, fixed now: if both measurements have
    equal, independent noise, their Pearson correlation is the share of the
    variation that is real, and so the highest R² any model could reach
    against one measurement. This is an estimate with assumptions, stated
    wherever it is used. Ipatasertib was measured once, so its own ceiling
    cannot be estimated; the two related drugs are a guide only.
