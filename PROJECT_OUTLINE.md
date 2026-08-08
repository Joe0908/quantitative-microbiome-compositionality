# Comprehensive project and training outline

## 1. Research question and evidential scope

The project asks:

> When the same published quantitative stool-microbiome profiles are analyzed
> as QMP, QMP-derived row-closed abundance, or CLR abundance, how stable are
> taxon-level inference and cross-validated disease discrimination?

The central claim is intentionally limited: abundance representation and CLR
zero handling can materially change feature-level findings, whereas QMP does
not show a consistent prediction advantage in these two cohorts. This is not
an equivalence or non-inferiority claim.

The study does not claim new CRC or IHD biomarkers, causal microbial effects,
same-disease external validation, or universal superiority of any
representation.

## 2. Datasets and primary contrasts

| Cohort | Contrast | Negative class | Positive class | N | Quantitative scale |
|---|---|---:|---:|---:|---|
| LCPM | CRC vs CTL | 205 CTL | 47 CRC | 252 | Species cells/g |
| MetaCardis | IHD372 vs MMC372 | 369 MMC | 303 IHD | 672 | Load-corrected MGS index |

LCPM also contains 337 adenoma participants used only in the source-aligned
three-group sensitivity analysis. MetaCardis contains overlapping published
analysis memberships; the primary IHD372 and MMC372 groups are participant-ID
disjoint.

The diseases and quantitative units differ. The two cohorts provide
cross-disease methodological replication, not biological replication, and raw
abundances or disease effects are never pooled.

## 3. Abundance representations

For quantitative matrix `X`, participant `i`, and taxon `j`:

1. **QMP** is the published quantitative value. LCPM supplies cells/g;
   MetaCardis supplies a load-corrected FPKM-derived index.
2. **QMP-derived row-closed abundance** is

   \[
   R_{ij}=\frac{X_{ij}}{\sum_k X_{ik}}.
   \]

   It is a deterministic mathematical transformation of QMP, not an
   independently measured native sequencing relative-abundance matrix. The
   QMP and row-closed matrices share upstream measurement error.
3. **Minimum-positive CLR** replaces a zero for taxon `j` with the minimum
   positive training-fold row-closed value for that taxon, takes natural logs,
   and subtracts the participant mean log abundance.
4. **Multiplicative CLR sensitivity** first closes the selected features,
   sets `delta = 1 / p²`, replaces every zero with `delta`, rescales non-zero
   entries by `1 - m_i × delta`, and applies CLR. Here `p` is the number of
   selected taxa and `m_i` is the number of zeros in participant `i`.

The minimum-positive rule is retained for continuity with the original
analysis. Multiplicative replacement tests whether conclusions depend on that
zero-handling choice.

## 4. Acquisition, construction, and QC

### Part 01 — acquire immutable inputs

- Download the two official supplementary workbooks through HTTPS.
- Store them under `data/raw/`.
- Compare each file with the fixed SHA-256 in `compositionality/config.py`.
- Stop before analysis if a checksum differs.

### Part 02 — construct LCPM

- Read participant metadata, the clean-species list, and species-level QMP.
- Preserve source participant and feature order.
- Verify 589 participants × 676 source features.
- Verify diagnosis counts: 205 CTL, 337 ADE, and 47 CRC.
- Retain 336 taxonomically clean candidate features.
- Derive row-closed abundance from the complete QMP row.
- Record total load, detected-feature count, and row-sum QC.

### Part 03 — construct and audit MetaCardis

- Stream the required worksheet XML to avoid loading the workbook's large
  formatting model.
- Audit 1,882 published membership rows representing 1,087 unique IDs.
- Confirm duplicate rows for an ID are identical apart from membership label.
- Create explicit HC, MMC, UMMC, IHD, ACS, CIHD, and HF membership indicators.
- Retain one quantitative profile per ID.
- Verify 994 complete profiles × 729 source MGS features.
- Map each MGS column to taxonomy and retain 416 clean bacterial candidates.
- Encode age, BMI, sex, nationality, diabetes, and four broad medication
  categories.
- Verify 303 IHD and 369 MMC profiles with no cross-group ID overlap.

## 5. Feature-level association analyses

### 5.1 LCPM primary analysis

Question: among CRC and CTL participants, do taxa differ under each abundance
representation?

1. Restrict to the 252 primary participants.
2. Calculate detection prevalence across all 252 participants without using
   CRC/CTL labels.
3. Retain taxa detected at row-closed abundance ≥ 1×10⁻⁶ in at least 5% of the
   pooled participants. This yields 93 taxa.
4. Analyze QMP, row-closed abundance, minimum-positive CLR, and multiplicative
   CLR separately.
5. For each taxon, run a two-sided Mann–Whitney U test.
6. Report rank-biserial effect size oriented as CRC minus CTL.
7. Apply Benjamini–Hochberg correction separately within each representation
   and zero-replacement specification.

These tests are unadjusted. Participant-level BMI, stool moisture, and
calprotectin needed for confounder adjustment are unavailable in the public
supplement.

### 5.2 LCPM source-aligned sensitivity

The published/source-aligned rule retains a taxon if it reaches 5% prevalence
in any of CTL, ADE, or CRC, yielding 112 taxa. Because diagnosis groups enter
feature selection, this result is retained as a sensitivity analysis rather
than the primary inferential result. Three-group 1%, 5%, and 10% analyses are
also retained for methodological context.

### 5.3 MetaCardis association models

Question: are IHD/MMC associations stable across prevalence, positive QMP,
positive row-closed abundance, and two CLR zero treatments?

1. Calculate pooled prevalence across all 672 participants without using the
   IHD/MMC label.
2. Retain 404 taxa present in at least 5% of participants.
3. For each taxon fit five components:

   - presence/absence logistic regression;
   - natural-log positive QMP OLS;
   - natural-log positive row-closed OLS;
   - minimum-positive CLR OLS across all participants;
   - multiplicative CLR OLS across all participants.

4. Positive-abundance models use only participants in whom the taxon is
   present.
5. OLS models use HC3 robust standard errors and finite-sample t inference.
6. Run three covariate specifications:

   - unadjusted disease indicator;
   - core: age, BMI, BMI-missing indicator, sex, nationality, and diabetes;
   - core plus antidiabetic, antihypertensive, lipid-lowering, and proton-pump
     inhibitor categories.

7. Remove constant or linearly dependent nuisance columns inside a sparse
   positive subset while retaining the disease term.
8. Apply BH correction separately within every component and covariate
   specification.

The medication-expanded specification is a sensitivity analysis under an
alternative covariate set. It is not a causal estimate of direct drug effects;
medication can be a confounder, mediator, or proxy for disease severity.

## 6. Prediction design

### 6.1 Fixed global settings

| Setting | Value |
|---|---|
| Splitter | Repeated stratified cross-validation |
| Folds × repeats | 5 × 10 |
| Random seed | 531 |
| Candidate taxa | 336 LCPM; 416 MetaCardis |
| Training filter | ≥5% in either training class |
| LCPM detection | Row-closed ≥1×10⁻⁶ |
| MetaCardis detection | QMP >0 |
| Primary estimator | `HistGradientBoostingClassifier` |
| Learning rate | 0.05 |
| Iterations | 150 |
| Maximum leaf nodes | 15 |
| Minimum samples per leaf | 10 |
| L2 regularization | 1.0 |
| Class handling | Balanced training sample weights |
| Early stopping | Disabled |
| Hyperparameter tuning | None |

The class-aware feature filter is a supervised prediction step and is valid
because it is fitted exclusively inside each training fold. It must not be
confused with the label-blind filter used for taxon-level hypothesis testing.

### 6.2 Exact operations within every fold

#### Step 1 — shared split

- Preserve published quantitative-matrix row order.
- Obtain train/test positions from the shared repeated-stratified splitter.
- Use the same positions for QMP, row-closed, both CLR variants, clinical, and
  combined models.
- Record participant and class counts.

#### Step 2 — training-fold taxon selection

- Start from the clean candidate set.
- Calculate detection prevalence in each class using training participants
  only.
- Retain a taxon if it reaches 5% in either training class.
- Record every selected taxon and selection frequency.

#### Step 3 — build abundance matrices

- Subset QMP and QMP-derived row-closed abundance to the same selected taxa.
- Keep zeros in QMP and row-closed models.
- Do not use test outcomes or test prevalence.

#### Step 4 — minimum-positive CLR

- For every selected taxon, learn its minimum positive row-closed value from
  the training fold.
- Replace training zeros and transform to CLR.
- Reuse the exact training replacement vector for held-out zeros.
- Record minimum and maximum fold replacement values.

#### Step 5 — multiplicative CLR sensitivity

- Use the same selected taxa.
- Apply the fixed `delta = 1 / p²` multiplicative replacement independently to
  each composition; no statistic is learned from held-out participants.
- Record `delta` for every fold.

#### Step 6 — fit microbiome models

- Fit separate fixed HGB models to QMP, row-closed, minimum-positive CLR, and
  multiplicative CLR.
- Use identical training labels and balanced weights.
- Save held-out positive-class probabilities.

#### Step 7 — construct MetaCardis clinical design

The clinical design contains age, BMI, BMI missingness, sex, France-vs-Denmark
nationality, diabetes, and four medication categories. For each fold:

- learn the BMI median from training participants only;
- impute training and held-out BMI with that value;
- retain the BMI-missing indicator.

#### Step 8 — fit estimator-matched clinical comparisons

- Fit clinical-only HGB using the same HGB settings and folds as the combined
  models.
- Fit QMP+clinical, row-closed+clinical, and minimum-positive CLR+clinical HGB.
- Retain a balanced logistic clinical-only model as contextual sensitivity,
  not as the estimator-matched incremental baseline.

#### Step 9 — fold and repeat metrics

For each held-out fold, save ROC AUC, average precision, balanced accuracy at
0.5, and Brier score. For each repeat, concatenate all five held-out folds so
every participant has exactly one out-of-fold prediction, then calculate:

- ROC AUC;
- average precision;
- balanced accuracy;
- sensitivity and specificity at 0.5;
- Brier score.

The 0.5 threshold is not optimized.

#### Step 10 — descriptive repeated-CV summary

- Report the mean, sample SD, literal minimum, and literal maximum across the
  ten repeat-level metrics.
- For paired models, report mean/SD/min/max of the ten within-repeat
  differences.
- Do not calculate a t-based interval, paired Wilcoxon P value, or BH q value.

Repeated-CV values are correlated because participants are reused. The
repeat-level range is descriptive resampling variability, not a confidence
interval. These results do not establish equivalence or non-inferiority.

## 7. Cross-cohort synthesis

Allowed operations:

- compare representations within each cohort;
- calculate within-cohort effect correlations, direction agreement, and call
  overlap;
- normalize strict species names and retain exactly one feature per species in
  each cohort;
- standardize effects within cohort and representation;
- compare representation-sensitivity gaps descriptively.

Prohibited operations:

- pooling CRC and IHD participants;
- pooling LCPM cells/g with the MetaCardis index;
- meta-analyzing different disease estimands;
- treating MetaCardis as biological validation of CRC;
- assigning independent-taxon P values to cross-species correlations.

After the revised pooled filters, exact one-to-one harmonization yields 49
shared species. Correlations across these taxa are descriptive because taxa
are not independent observations.

## 8. Outputs and audit trail

| Output group | Content |
|---|---|
| Processed data | Metadata, taxonomy, QMP, row-closed abundance, QC JSON |
| Associations | Effects, SEs where applicable, P values, BH q values, fit status |
| CV repeats | One row per cohort/model/repeat with complete OOF metrics |
| CV folds | Held-out metrics for every model and fold |
| Fold audit | Counts, selected taxa, CLR replacement range/delta, BMI median |
| Selection frequency | Number and fraction of 50 folds selecting each candidate |
| Sensitivity | LCPM filter comparison and CLR zero-replacement comparison |
| Synthesis | Descriptive CV differences, call counts, effect stability, 49 species |

## 9. Acceptance criteria

A run is accepted only when:

- source dimensions and diagnosis/membership counts match the fixed checks;
- all row-closed profiles sum to one within numerical tolerance;
- no association filter uses disease labels in the primary analysis;
- all prediction preprocessing is learned within the training fold;
- all held-out probabilities are finite and complete;
- repeated-CV exports contain no inferential CI, P, or q fields;
- association counts, AUC checkpoints, effect correlations, and shared-species
  count match `expected_results.json`;
- unit tests pass; and
- no output pools incompatible diseases, abundance scales, or estimands.

## 10. Final interpretation boundary

The revised evidence supports a narrower and more defensible conclusion than
the first analysis. LCPM feature discoveries were highly sensitive to whether
the tested taxa were selected with or without diagnosis labels. CLR results
and prediction also changed with zero replacement. Across both cohorts, QMP
did not consistently improve disease discrimination, and MetaCardis clinical
variables dominated microbiome-only prediction. Quantitative measurement is
valuable, but it does not by itself eliminate compositional sensitivity,
outcome-dependent filtering bias, or clinical confounding.
