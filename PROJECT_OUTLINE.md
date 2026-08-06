# Comprehensive project outline

## 1. Project question and final scope

The project asks one focused methodological question:

> When the same stool microbiome profiles are represented as quantitative
> abundance (QMP), row-closed relative abundance (RMP), or centered log-ratio
> abundance (CLR), how much do disease prediction and taxon-level conclusions
> change?

The scope is deliberately limited to two suitable published datasets:

1. **Galazzo/LCPM** supplies a true species-level cells/g benchmark for
   colorectal cancer.
2. **MetaCardis** supplies an independent cardiovascular cohort with a
   microbial-load-corrected quantitative abundance index and richer metadata.

This is a data-replacement and robustness project, not a new sequencing
pipeline. The published quantitative matrices are processed into analysis-ready
datasets; synthetic data are not generated.

## 2. Prespecified analyses

### 2.1 Primary disease contrasts

| Cohort | Contrast | Negative class | Positive class | N |
|---|---|---:|---:|---:|
| Galazzo/LCPM | CRC vs CTL | 205 CTL | 47 CRC | 252 |
| MetaCardis | IHD372 vs MMC372 | 369 MMC | 303 IHD | 672 |

MetaCardis published labels overlap in several other designs. The IHD372 and
MMC372 groups are participant-ID disjoint, so they are used as the primary
contrast. Nested ACS, CIHD, and HF labels are not added as extra cases.

### 2.2 Representations

For a participant-by-taxon quantitative matrix \(X\):

1. **QMP:** the published quantitative value, with structural zeros retained.
2. **RMP:** row closure,
   \[
   R_{ij}=\frac{X_{ij}}{\sum_k X_{ik}}.
   \]
3. **CLR:** after feature selection, replace a zero for feature \(j\) with the
   minimum positive training-fold value for that feature, take natural logs,
   and subtract the participant's mean log abundance,
   \[
   \operatorname{clr}(R_{ij})=\log(R_{ij}^{*})-
   \frac{1}{p}\sum_{k=1}^{p}\log(R_{ik}^{*}).
   \]

For held-out data, the replacement values learned from the training fold are
reused unchanged.

## 3. Data-processing phases

### Phase 01 — acquisition and provenance

**Input:** two public Nature Medicine supplementary workbooks.

**Actions:**

- Download through HTTPS.
- Cache under `data/raw/`.
- Verify the fixed SHA-256 before analysis.
- Stop on any checksum mismatch.

**Output:** verified immutable input workbooks.

### Phase 02 — Galazzo/LCPM construction

**Source sheets:** participant metadata, clean-species list, and species QMP
cells/g matrix.

**Actions:**

- Standardize participant and diagnosis fields.
- Preserve the published participant and feature order.
- Validate 589 participants and 676 source taxa.
- Retain the prespecified 336 taxonomically clean candidate taxa.
- Derive RMP by closing each full QMP row.
- Add sample-level load and detected-feature audit fields.
- Confirm diagnosis counts: 205 CTL, 337 ADE, and 47 CRC.

**Outputs:** metadata, taxonomy, QMP cells/g, RMP, and QC JSON.

### Phase 03 — MetaCardis construction and overlap audit

**Source sheets:** phenotype/medication tables, quantitative MGS matrix, and
taxonomy.

**Actions:**

- Stream the required worksheet XML to avoid loading the workbook's large style
  model into memory.
- Audit 1,882 published label rows representing 1,087 unique participant IDs.
- Verify that duplicated rows for an ID differ only in the published `Status`
  membership.
- Create explicit Boolean membership fields for HC275, MMC372, UMMC, IHD372,
  ACS, CIHD, and HF.
- Keep one quantitative profile per ID.
- Validate 994 complete QMP/RMP profiles and 729 source MGS features.
- Map MGS columns one-to-one to taxonomy and retain 416 clean bacterial
  candidates.
- Convert sex, diabetes, and four broad drug categories to numeric fields.
- Confirm 303 IHD and 369 MMC profiles with no cross-group ID overlap.

**Outputs:** audited metadata, taxonomy, QMP index, RMP, and QC JSON.

## 4. Model-training specification

### 4.1 Global fixed settings

| Setting | Value |
|---|---|
| Splitter | Repeated stratified cross-validation |
| Folds × repeats | 5 × 10 |
| Random seed | 531 |
| Primary prevalence | 5% in either training class |
| LCPM detection threshold | RMP ≥ 1×10⁻⁶ |
| MetaCardis detection threshold | QMP > 0 |
| Primary classifier | `HistGradientBoostingClassifier` |
| Learning rate | 0.05 |
| Iterations | 150 |
| Maximum leaf nodes | 15 |
| Minimum samples per leaf | 10 |
| L2 regularization | 1.0 |
| Class handling | Balanced training weights |
| Early stopping | Disabled |

There is no hyperparameter search. Fixing the model before comparison prevents
representation-specific tuning from becoming an alternative explanation for a
QMP/RMP/CLR difference.

### 4.2 Exact content of every training fold

The following steps are repeated 50 times per cohort.

#### Training step 1 — define the fold

- Read the train/test positions from the shared repeated-stratified splitter.
- Preserve the published quantitative-matrix row order.
- Use exactly the same positions for every representation and optional
  baseline.
- Record train/test counts and class counts.

#### Training step 2 — fit the feature filter

- Start from 336 LCPM or 416 MetaCardis clean candidates.
- Using **training participants only**, calculate detection prevalence within
  each binary class.
- Retain a taxon if it reaches 5% prevalence in either class.
- Record the number and identity of selected taxa.

This produces a varying feature set by fold. At the reference run, the mean
number selected is 116.58 for LCPM and 409.56 for MetaCardis.

#### Training step 3 — construct fold matrices

- Subset QMP and the already row-closed RMP to the selected taxa.
- RMP closure itself is a deterministic within-participant calculation and uses
  no outcome or other participant.
- Keep zeros in QMP and RMP; do not add a prediction pseudocount.

#### Training step 4 — fit the CLR replacement reference

- For each selected taxon, find its minimum positive RMP in the training fold.
- Store that vector as the fold-specific replacement reference.
- Replace training zeros, log, and center within participant.
- Apply the **same stored vector** to held-out zeros before held-out log-centering.
- Record minimum and maximum replacement values for leakage auditing.

#### Training step 5 — calculate class weights

- For a training fold with \(n\) participants and \(n_c\) participants in class
  \(c\), assign weight \(n/(2n_c)\).
- Fit all three microbiome models with the same training labels and weights.

#### Training step 6 — fit the QMP model

- Input: selected raw quantitative features.
- Estimator: fixed histogram-gradient boosting model.
- Output: held-out probability of CRC for LCPM or IHD for MetaCardis.

#### Training step 7 — fit the RMP model

- Input: the same selected taxa expressed as row-closed relative abundance.
- Estimator and weights: identical to the QMP model.
- Output: held-out positive-class probabilities.

#### Training step 8 — fit the CLR model

- Input: fold-local CLR matrix for the same selected taxa.
- Estimator and weights: identical to QMP and RMP.
- Output: held-out positive-class probabilities.

#### Training step 9 — optional MetaCardis clinical baseline

The clinical design contains:

- age;
- BMI and a BMI-missing indicator;
- male indicator;
- France versus Denmark nationality indicator;
- diabetes;
- antidiabetic medication;
- antihypertensive medication;
- lipid-lowering medication; and
- proton-pump inhibitor use.

Eight BMI values are missing. Within every fold, the training-fold median is
learned and applied to both training and held-out BMI while the missingness
indicator is retained. The full primary-cohort audit median is 27.568 kg/m².
A balanced L2 logistic regression produces the held-out clinical probability.

#### Training step 10 — optional combined models

- Append the ten clinical fields to QMP, RMP, or CLR.
- Fit three additional fixed histogram-gradient boosting models.
- These models are secondary context checks; they are not used to decide
  whether QMP solves compositionality.

#### Training step 11 — score the held-out fold

For every fitted model, save:

- ROC AUC;
- average precision;
- balanced accuracy at probability 0.5; and
- Brier score.

The 0.5 threshold is not optimized.

#### Training step 12 — assemble repeat-level out-of-fold predictions

- Concatenate the five held-out folds so every participant has one prediction
  for that repeat.
- Calculate ROC AUC, average precision, balanced accuracy, sensitivity,
  specificity, and Brier score on the full out-of-fold vector.
- Repeat for all ten repeats.

#### Training step 13 — summarize and compare models

- Report mean, sample SD, and a t-based exploratory 95% interval across the ten
  repeat metrics.
- Compare representation pairs using a paired two-sided Wilcoxon test on the
  ten repeat values.
- Apply BH correction within the comparison family.

These intervals and tests are descriptive because repeated-CV estimates are
correlated. They are not substitutes for external validation.

## 5. Differential-association analyses

### 5.1 Galazzo/LCPM

- Define the 5% analysis feature set across CTL, ADE, and CRC, then test the
  prespecified CRC-vs-CTL contrast.
- Use two-sided Mann–Whitney U tests for QMP, RMP, and CLR.
- Report rank-biserial effect size oriented as CRC minus CTL.
- Apply BH correction separately to each representation.
- Repeat global CTL/ADE/CRC sensitivity analyses at 1%, 5%, and 10% prevalence.

Primary fixed set: 112 taxa. Expected calls: 8 QMP, 8 RMP, and 1 CLR.

### 5.2 MetaCardis

Use 410 fixed 5%-prevalence features and separate presence from positive
abundance:

1. **Prevalence component:** logistic regression of presence on IHD status.
2. **QMP non-zero component:** OLS of natural-log positive QMP with HC3 robust
   standard errors and finite-sample t tests.
3. **RMP non-zero component:** the same model on natural-log positive RMP.
4. **CLR component:** OLS on CLR across all participants with HC3 robust
   standard errors and finite-sample t tests.

Specifications are:

- unadjusted disease indicator;
- core adjusted: age, BMI, BMI missingness, sex, nationality, diabetes; and
- core plus the four medication categories.

Constant or collinear nuisance indicators are removed inside a taxon's positive
subset without dropping the disease term. BH correction is applied separately
within every component and adjustment specification.

## 6. Cross-cohort synthesis

### What is allowed

- Compare QMP/RMP/CLR performance **within** each cohort.
- Compare representation-specific effect correlation and direction agreement
  within each cohort.
- Normalize strict two-token species names and require exact matches.
- Apply the completed project's conservative MetaCardis cross-cohort flag:
  classified bacterial labels containing neither `sp.` nor the substring
  `bacterium`; then apply the strict two-token binomial rule.
- Require exactly one feature per species in each cohort.
- Standardize effects within cohort and representation across the full tested
  feature set.
- Compare QMP-minus-RMP and QMP-minus-CLR standardized gaps for the 51 shared
  species.

### What is prohibited

- Pooling CRC and IHD participants.
- Pooling LCPM cells/g with the MetaCardis abundance index.
- Meta-analysis of disease effects with different estimands and adjustments.
- Interpreting shared-species gap correlation as shared disease biology.
- Arbitrarily choosing or summing one of several MetaCardis MGS features with
  the same species name; duplicated labels are excluded.

## 7. Output and audit plan

| Output group | Specific content |
|---|---|
| Processed LCPM | metadata, taxonomy, QMP, RMP, QC checks |
| Processed MetaCardis | overlap-aware metadata, taxonomy, QMP, RMP, QC checks |
| Associations | full effects, SEs where applicable, p values, BH q values, fit status |
| CV repeats | one row per cohort/model/repeat with full OOF metrics |
| CV folds | held-out metrics for every fold and model |
| Fold audit | sample counts, selected-feature count, CLR replacement range |
| Selection frequency | number and fraction of 50 folds selecting each candidate |
| Synthesis | headline JSON, CV comparison, DA call counts, effect stability |
| Shared species | one-to-one species effects, z scores, and representation gaps |

## 8. Acceptance criteria

The run is accepted when all of the following hold:

- LCPM is 589 × 676 before primary subsetting, with 336 clean candidates.
- MetaCardis is 994 × 729 before primary subsetting, with 416 clean candidates.
- Primary groups are 47/205 and 303/369 and contain no MetaCardis cross-group ID
  overlap.
- Every RMP row sums to one within numerical tolerance.
- Training-fold filtering and CLR replacement are audited for all 100 cohort
  folds.
- Primary AUCs match `expected_results.json` within the documented software
  tolerance.
- Differential call counts match the expected checkpoints.
- Exact harmonization yields 51 shared one-to-one species.
- No output pools incompatible disease effects or quantitative scales.

## 9. Final interpretation

The replacement data solve the original design limitation: the same samples can
now be examined as quantitative, relative, and log-ratio profiles. The results
do not support a universal prediction advantage for QMP, but they show that
taxon discovery and uncertainty can depend materially on representation.
MetaCardis also shows that clinical and medication structure can matter more
than the small QMP-versus-RMP AUC difference. Quantitative measurement is
therefore valuable, but it does not remove the need for compositional
sensitivity analysis or careful confounder data.
