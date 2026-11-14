# Drug response beyond a single marker: report

**Author.** Paul McNeill. Written on 13 November 2026 from the logged results in this repository; not yet peer reviewed. Every number below comes from a table in `results/`.

**How this work was produced.** Both studies were built making the use of AI tools, under research rules I set; each study’s log records how. I re-ran both on my own hardware, reviewed every logged decision, and take responsibility for the work and its conclusions.

## 1. Question

For a drug where the status of one gene (PTEN) is a known but incomplete
guide to response, how much more of the response does wider molecular data
explain, and are the patterns found stable enough to trust?

The setting mirrors, on public data, a common problem in treatment-response
research: one marker that does not explain everything, several kinds of
measurement that do not cover every sample, and subgroups that may be
unstable or technical. No prediction was recorded before the results existed.

## 2. Data and setup

- **Data.** 924 cancer cell lines from the GDSC2 drug screen, with mutations,
  copy number, gene expression and protein abundance from the Sanger Cell
  Model Passports. 128 cell lines have PTEN altered.
- **Drug.** Ipatasertib, which blocks AKT, the protein that PTEN loss
  switches on. Three related drugs are used for replication.
- **Response.** AUC, the area under the dose-response curve.
- **Background factors.** Cancer type and sequencing centre, because either
  could pass for signal.
- **Validation.** 5-fold cross-validation repeated 5 times on saved folds.
  Every step that learns from data happens inside the training folds.
- **Protocol.** Written before any model was fitted. Ten dated amendments,
  four of them made after seeing results and labelled as such.

## 3. Methods

1. **A ladder of baselines**, each adding one kind of information: PTEN
   status; cancer type and centre; both; then 2,000 genes with an elastic net
   and with gradient boosting.
2. **A leak, measured.** Genes chosen by correlation with response using all
   cell lines (wrong) or inside each fold (right), on real and on shuffled
   response.
3. **A second data type with gaps.** Protein data hidden at random for some
   cell lines; four ways of handling that.
4. **Cluster stability** among PTEN-altered cell lines, against reference
   data with no groups, before and after removing tissue and centre.
5. **Explanation stability.** How far the top-20 feature lists of refitted
   models overlap.
6. **Repeatability of the response.** Agreement between two measurements of
   the same drug on the same cell lines.

## 4. Results

**What each kind of information adds (held-out R²).**

| Inputs | R² |
| --- | ---: |
| PTEN status | 0.051 |
| Cancer type + sequencing centre | 0.068 |
| PTEN + background | 0.112 |
| + 2,000 genes, elastic net | 0.144 |
| + 2,000 genes, gradient boosting | 0.187 |

- PTEN is a real marker and a weak one. It still adds about 4 points on top
  of cancer type.
- Expression adds +0.033 with the elastic net (95% interval +0.005 to
  +0.060) and +0.076 with boosting (+0.014 to +0.130).
- Boosting scores higher than the elastic net, but the difference (+0.043,
  95% interval -0.003 to +0.086) is not clear of zero.
- More than 80% of the response is unexplained by every model.

**The leak.** Choosing genes with all cell lines gives R² 0.210 on the real
response and 0.137 on a shuffled response that contains nothing. Done inside
the folds, the same recipe gives -0.005 and -0.173. The false 0.137 is as
large as the real result of the elastic net.

**Protein data added nothing.** Late fusion matched expression alone (0.152
to 0.161 against 0.156) at every level of overlap. Training only on complete
cases, or filling gaps with medians, was worse, and collapsed when 25% of
cell lines had protein.

**No stable subgroups beyond tissue were found with these methods.** Three clusters are perfectly stable
as published (stability 1.00 against 0.68 for data with no groups), and 94%
of cell lines sit in the commonest cluster of their own tissue. After
removing tissue and centre, no clustering is both stable and clearly above
the reference. Data with no groups scores 0.80 to 0.87 when cut in two, so
stability alone proves little.

**Explanations are not stable.** Top-20 lists from 25 refits overlap by 0.20
on average. Only PTEN status and one gene are in the top 20 at least 80% of
the time.

**Replication is partial.** The gain from expression is clear for one of the
three related drugs.

**Noise is not the main reason for the gap** (exploratory). A related drug
measured twice at the same doses agrees at 0.81, so about 80% of its
variation is repeatable. The best model for that drug explains 0.18.

## 5. Against the criteria fixed in advance

| Criterion for a negative result | Outcome |
| --- | --- |
| Expression gain interval includes zero | Not met for Ipatasertib. Met for 2 of the 3 related drugs |
| No clustering reaches stability 0.50 after adjustment | Not met as written, but the criterion had no reference. Against the reference, nothing qualifies |
| Top-20 overlap below 0.30 | Met: 0.20 |

## 6. What failed, and what is uncertain

- Gains were first computed on predictions averaged over repeats, which
  flattered the less stable model (+0.094 became +0.076). Corrected.
- Feature importance was first averaged only over the fits where a gene was
  an input. Corrected.
- A sanity check declared in advance was nearly empty: 10 of the 12 pathway
  genes it named were filtered out before the models saw them.
- The honestly done version of the leak recipe scored about zero on the real
  response, because a few extreme predictions wrecked R². This was not
  expected.

## 7. Limitations

- Cell lines are not patients. Nothing here is about treatment.
- One screen from one laboratory. The related drugs share its cell lines.
- The intervals resample cell lines only, so real uncertainty is larger.
- Two thirds of cell lines barely respond, so R² is decided by a minority.
- Missing protein data was simulated at random. Models were not tuned.
- The noise ceiling comes from a related drug, not from Ipatasertib.
- No biological interpretation has been checked by a biologist.

## 8. Related work

Not written. It needs the papers to be read in full first. `CITATION.md`
lists the dataset papers; only their journal details were checked.

## 9. What would come next

- An independent dataset for the same drugs.
- Leave-one-tissue-out validation, which asks a harder question.
- Stability selection, to get a gene list that survives refitting.
- A knowledge-based gene set beside the variance filter.
- A learned model of missing protein data in place of medians.
