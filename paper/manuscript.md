# A paired benchmark of quantitative, row-closed, and log-ratio gut microbiome profiles in two public cohorts

**Short title:** Quantitative versus compositional microbiome representations

**Authors:** [AUTHOR NAME(S) TO CONFIRM]

**Affiliations:** [AFFILIATION(S) TO CONFIRM]

**Corresponding author:** [NAME, POSTAL ADDRESS, AND EMAIL TO CONFIRM]

## Abstract

### Background

Relative-abundance microbiome profiles are compositional, whereas quantitative microbiome profiling (QMP) retains microbial-load information. Whether that information consistently changes feature-level associations or disease discrimination in public human cohorts remains uncertain.

### Methods

We performed a paired reanalysis of two public gut-microbiome datasets with microbial-load information. The primary LCPM contrast comprised 205 lesion-free controls and 47 participants with colorectal cancer (CRC); the primary MetaCardis contrast comprised 369 metabolically matched controls (MMC) and 303 participants with ischaemic heart disease (IHD). Within each cohort, published quantitative values were compared with a deterministically derived row-closed composition and centered log-ratio (CLR) representations using minimum-positive and multiplicative zero replacement. Association filtering was based on pooled prevalence without outcome labels; the source-aligned LCPM group-union rule was retained as a sensitivity analysis. Disease discrimination used fixed histogram-gradient-boosting models under 10 repeats of stratified five-fold cross-validation, with identical splits and training-fold-only preprocessing. Repeated-cross-validation summaries were treated as descriptive because repeats reuse participants.

### Results

Quantitative representation did not show a consistent discrimination advantage. Mean ROC AUCs for QMP, row-closed abundance, minimum-positive CLR, and multiplicative-replacement CLR were 0.659, 0.653, 0.642, and 0.663 in LCPM and 0.639, 0.647, 0.643, and 0.652 in MetaCardis. The descriptive mean QMP-minus-row-closed difference was +0.006 in LCPM (repeat range, −0.053 to +0.065) and −0.008 in MetaCardis (−0.032 to +0.021). Under pooled outcome-blind filtering, LCPM tested 93 species and identified one feature under each representation; the source-aligned group-union rule tested 112 species and produced 8 QMP, 8 row-closed, and one call under each CLR rule. Core-adjusted MetaCardis models tested 404 features and identified 10 prevalence, 6 QMP, 0 row-closed, 13 minimum-positive CLR, and 14 multiplicative-CLR features. With medication covariates, the corresponding counts were 0, 0, 1, 3, and 0. An estimator-matched MetaCardis clinical-only model achieved an AUC of 0.894, compared with 0.896 for row-closed abundance plus clinical variables.

### Conclusions

Across two disease settings, QMP did not consistently improve internally cross-validated discrimination. Feature-level conclusions depended on outcome-blind filtering, abundance representation, CLR zero replacement, and covariate specification. These results support paired sensitivity analyses rather than a universal hierarchy of representations. The cross-disease comparison is methodological replication, not biological validation or evidence of equivalence.

**Keywords:** quantitative microbiome profiling; compositional data; centered log-ratio; colorectal cancer; ischaemic heart disease; machine learning; differential abundance; microbial load

## Introduction

Sequencing-based microbiome studies usually report each taxon as a fraction of the reads observed in a sample. These profiles are constrained to a constant total and are therefore compositional: the measured abundance of one taxon depends mathematically on all other measured taxa [1]. This dependence complicates interpretation because relative changes may reflect altered total microbial load, altered abundance of another organism, or both. Consequently, disease-associated relative-abundance signals cannot automatically be interpreted as changes in the number of microbial cells.

Quantitative microbiome profiling combines community composition with an independent estimate of microbial load. Flow-cytometry-based QMP has shown that microbial load can explain biologically important variation that is obscured by relative profiles [2]. Experimental benchmarking has likewise indicated that quantitative approaches can reduce compositional and sampling-depth biases in association analyses [3]. Quantitative values, however, are not methodologically uniform: cells per gram, gene-count-normalized indices, spike-in estimates, and other load-corrected measurements have distinct units and error structures. Quantification therefore adds information but does not remove the need for careful normalization, confounder control, and sensitivity analysis.

The choice of abundance representation is especially consequential for feature-level inference. Differential-abundance methods can produce markedly different discoveries on the same microbiome data [4], and recent realistic benchmarks have emphasized the joint influence of method choice and confounding [5]. Prediction poses a related but distinct question. A transformation may change individual coefficients and false-discovery-rate calls without materially changing multivariable discrimination, particularly when the classifier can exploit correlated features. Machine-learning performance can also be inflated by data leakage, outcome-informed preprocessing, or unpaired validation designs [8]. Recent disease-classification benchmarks further suggest that the interaction between normalization, feature selection, and classifier choice is context dependent [9].

Two published cohorts provide an opportunity for a real-data, paired comparison. The LCPM colorectal-neoplasia study released species-level quantitative values in cells per gram for 589 participants and reported strong effects of microbial load and confounding on candidate CRC markers [6]. MetaCardis released load-corrected metagenomic species (MGS) profiles and detailed clinical and medication metadata across the cardiometabolic disease spectrum [7]. These cohorts differ in disease, measurement scale, and covariate availability, so their raw abundances and disease effects cannot be pooled. They can, however, independently test whether the same abundance-representation choice has similar methodological consequences.

We therefore asked: when the same participants and microbial features are represented as published quantitative abundance, a QMP-derived row-closed composition, or CLR abundance, how much do disease discrimination and feature-level conclusions change? We expected that QMP would not provide a universal predictive advantage and assessed, without presuming direction, how filtering, zero replacement, and covariate specification affected feature-level conclusions. The design used paired, leakage-controlled cross-validation; outcome-blind association filtering; two CLR zero-replacement rules; medication sensitivity analysis; an estimator-matched clinical comparison; a guarded shared-species analysis; and a post hoc analysis of total microbial load.

## Methods

### Study design and reporting framework

This was a secondary analysis of anonymized, publicly available supplementary data from two previously published human gut-microbiome studies [6,7]. No new participants were recruited and no new specimens were collected. The original publications describe participant consent and institutional ethical approvals; the present analysis used only public, de-identified tables. Reporting was organized with reference to the STORMS checklist for human microbiome research [12].

All analyses were performed separately within cohort. LCPM and MetaCardis were not concatenated because their diseases, study designs, taxonomic features, and quantitative-abundance units differ. Cross-cohort synthesis was restricted to methodological comparisons and standardized within-cohort representation gaps.

### LCPM data source and participants

LCPM data were obtained from Supplementary Tables S1, S6, and S14 of Tito et al. [6]. The public supplement contained 589 participants: 205 lesion-free controls (CTL), 337 participants with adenoma or colorectal polyps (ADE), and 47 participants with CRC. Supplementary Table S14 contained 676 species-level QMP features expressed as cells per gram. Sample identifiers aligned one-to-one between metadata and the quantitative matrix, and no negative QMP values were present.

The primary binary contrast was CRC versus CTL (n=252; 47 CRC and 205 CTL). ADE samples were retained only for the source-aligned three-group-union filtering sensitivity and prevalence-threshold analyses. The clean candidate list was derived from Supplementary Table S6 and comprised 336 species labels. Public participant-level BMI, faecal moisture, and calprotectin values were unavailable in the released supplement and therefore could not be included; all LCPM disease-association results are consequently unadjusted.

### MetaCardis data source, participant audit, and primary contrast

MetaCardis data were obtained from Supplementary Tables ST5, ST9, ST10, and ST14 of Fromentin et al. [7]. The public phenotype tables contained 1,882 cohort-label rows representing 1,087 unique participant identifiers. Repeated rows for a participant were identical except for published analytical-group membership. We therefore retained one quantitative profile per identifier and represented group memberships as explicit Boolean fields, avoiding treatment of repeated labels as independent samples.

Among 994 participants with complete quantitative profiles across 729 MGS features, the primary contrast used the participant-disjoint IHD372 and MMC372 groups: 303 participants with IHD and 369 MMC participants (n=672). Nested acute coronary syndrome, chronic IHD, and heart-failure labels were not added as additional independent cases. Taxonomy from ST5 was mapped one-to-one to matrix columns. Analyses began with 416 bacterial features that were not labelled as unclassified. Because multiple MGS features can map to the same species label, MetaCardis results are reported as feature-level rather than unique-species discoveries.

The MetaCardis quantitative matrix is a microbial-load-corrected FPKM-derived index, not taxon-specific cells per gram [7]. Its magnitude was therefore interpreted only within MetaCardis and was not compared numerically with LCPM QMP.

### Abundance representations

For each cohort, the published quantitative participant-by-feature matrix was denoted QMP, (X). A row-closed composition was derived for participant (i) and feature (j) as

\[
Rᵢⱼ = Xᵢⱼ / Σₖ Xᵢₖ.
\]

This matrix is referred to throughout as **row-closed abundance**. It is derived from the public QMP matrix and is not the original native sequencing relative-microbiome profile. This distinction prevents the derived composition from being misinterpreted as an independently measured RMP.

CLR values were calculated after feature selection under two zero-replacement rules. In the minimum-positive rule, each zero was replaced by that feature's minimum positive value in the relevant reference data. The transformed value was

\[
clr(R*ᵢⱼ) = ln(R*ᵢⱼ) − (1/p) Σₖ₌₁ᵖ ln(R*ᵢₖ),
\]

where (p) is the number of retained features. For association analyses, the minimum-positive reference was fitted on the complete analysis set. For prediction, feature-specific replacement values were fitted using the training fold only and applied unchanged to the held-out fold.

As a zero-handling sensitivity analysis, selected features were first re-closed to sum to one. For a row with (m) zeros, each zero was replaced by δ=1/p² and each non-zero component was multiplied by (1−mδ), preserving a positive unit-sum composition before CLR transformation [15]. This multiplicative rule uses only the row being transformed; prediction nevertheless retained fold-local feature selection.

### Feature filtering

In LCPM, detection was defined as row-closed abundance ≥1×10⁻⁶. The primary association filter retained a feature when detected in at least 5% of the pooled CRC-versus-CTL analysis set, without using outcome labels, yielding 93 tested species [14]. A source-aligned sensitivity retained features detected in at least 5% of any of CTL, ADE, or CRC, yielding 112 species; because that rule uses diagnostic-group membership, it was not treated as the confirmatory filter. Three-group sensitivity analyses repeated the source-aligned rule at 1%, 5%, and 10% prevalence. In each prediction fold, filtering used the pooled training set only, without training labels, and retained features with at least 5% prevalence.

In MetaCardis, detection was defined as QMP >0. The primary association filter retained features present in at least 5% of the pooled 672-participant analysis set, yielding 404 features. Prediction refitted the same pooled 5% rule within each training fold. Neither outcome labels nor held-out participants influenced feature inclusion.

### LCPM feature-level analyses

For each of the 93 primary features, CRC and CTL distributions were compared separately under QMP, row-closed, minimum-positive CLR, and multiplicative-replacement CLR representations using a two-sided asymptotic Mann–Whitney U test. The same analyses were repeated for the 112-feature source-aligned set. Rank-biserial correlation was reported as the effect size, oriented so that positive values indicate higher ranks in CRC. Because sparse features frequently had a median of zero in both groups, inference was based on the full rank distribution rather than median differences alone.

For the source-aligned CTL–ADE–CRC sensitivity analysis, QMP and row-closed values were tested with Kruskal–Wallis tests. CLR values were assessed using both one-way analysis of variance and Kruskal–Wallis tests to expose dependence on the parametric or rank-based specification. Benjamini–Hochberg correction [10] was applied within each filter specification, representation, prevalence threshold, and test family.

### MetaCardis two-part feature-level analyses

MetaCardis feature-level analysis separated prevalence from abundance among non-zero observations. For each of 404 features, presence/absence was modeled with binomial logistic regression. Among participants with non-zero abundance for that feature, log-QMP and log-row-closed values were analyzed with ordinary least squares using HC3 heteroskedasticity-robust standard errors and finite-sample t inference. Minimum-positive and multiplicative-replacement CLR values were analyzed across all participants with the same HC3 specification. Disease effects were oriented as IHD minus MMC.

Three covariate specifications were fitted. The unadjusted model included disease status only. The core-adjusted model additionally included age, BMI, a BMI-missingness indicator, sex, nationality, and diabetes. Eight missing BMI values were replaced by the median (27.568 kg/m²) while retaining the missingness indicator. The medication-sensitivity model added antidiabetic, antihypertensive, lipid-lowering, and proton-pump-inhibitor medication indicators. Medication may be a confounder, mediator, treatment-indication proxy, or marker of disease severity; this expanded model was therefore interpreted only as sensitivity to an alternative covariate set, not as a causal estimate of medication or a direct disease effect. Rank-deficient or invariant nuisance columns within a feature-specific positive subset were omitted explicitly. Benjamini–Hochberg correction was applied separately within each covariate specification and model component.

These models implement a transparent two-part strategy and should not be described as results from a named two-part software framework.

### Disease-discrimination models

Disease discrimination was evaluated separately for LCPM CRC versus CTL and MetaCardis IHD versus MMC. The primary classifier was scikit-learn's `HistGradientBoostingClassifier` [11] with fixed parameters: learning rate 0.05, 150 iterations, maximum 15 leaf nodes, minimum 10 samples per leaf, L2 regularization 1.0, early stopping disabled, and random seed 531. No hyperparameter search was performed. Training samples received balanced class weights.

We used repeated stratified five-fold cross-validation with 10 repeats and seed 531. QMP, row-closed, minimum-positive CLR, and multiplicative-replacement CLR models used identical participant splits. Within each fold, pooled prevalence filtering was fitted only on training samples; the same retained feature identities were used for all four representations. Minimum-positive CLR replacement values were fitted only on the training fold. Multiplicative replacement was calculated independently within each row after the fold-local feature set was fixed. The fitted models produced held-out probabilities for every participant once per repeat.

Primary model performance was summarized by ROC AUC. Average precision, balanced accuracy, sensitivity, specificity, and Brier score were retained as secondary metrics. Classification metrics used a fixed probability threshold of 0.5; no threshold was selected from the data. Means, sample standard deviations, and actual minimum-to-maximum ranges were calculated across the 10 repeat-level out-of-fold estimates. Because the repeats reuse participants and overlap in training data, they are not independent replicates. We therefore performed no t-based confidence interval, Wilcoxon test, or multiplicity correction across repeat-level metrics. These summaries neither establish statistical equivalence nor constitute a non-inferiority analysis.

### MetaCardis clinical and combined models

For MetaCardis, secondary clinical baselines used age, BMI, BMI missingness, sex, nationality, diabetes, and the four medication categories. BMI median imputation was fitted within each training fold. A balanced logistic-regression clinical baseline was retained for context. For the estimator-matched comparison, clinical-only, QMP-plus-clinical, row-closed-plus-clinical, and minimum-positive-CLR-plus-clinical models all used the same fixed histogram-gradient-boosting estimator and folds.

Estimator-matched repeat-level differences were summarized descriptively, with no inferential test for the same dependence reasons described above. No calibration curve, decision-curve analysis, externally fixed clinical threshold, or external validation was performed; these models do not establish clinical utility.

### Total microbial-load sensitivity analysis

As a post hoc exploratory analysis, total microbial load alone was compared between the two primary groups in each cohort. LCPM load was the row sum of the published cells-per-gram species matrix. MetaCardis load was the public microbial-load variable and was available for 298 IHD and 366 MMC participants. Groups were compared with two-sided Mann–Whitney U tests and rank-biserial effects. Discrimination was summarized as the ROC AUC obtained by using load directly as the score. Stratified percentile-bootstrap 95% confidence intervals used 10,000 class-wise resamples and seed 531.

### Cross-cohort synthesis

Exact one-to-one species harmonization was restricted to strict two-token binomial species labels that occurred once in each primary, outcome-blind feature set, yielding 49 shared species. Within each cohort and representation, feature effects were standardized to mean zero and unit variance. For each species, QMP-minus-row-closed and QMP-minus-minimum-positive-CLR standardized representation gaps were calculated. Pearson and Spearman correlations described the correspondence of these gaps across cohorts. No feature-level correlation P values were used because taxa are not independent replicates. These analyses assessed representation sensitivity; they did not pool CRC and IHD disease effects or claim biological replication.

### Software, code, and reproducibility

The reference analysis used Python 3.12.13, NumPy 2.3.5, pandas 2.2.3, SciPy 1.17.0, scikit-learn 1.8.0, statsmodels 0.14.6, and openpyxl 3.1.5. The repository contains supported version ranges, an exact `requirements-lock.txt`, fixed input checksums, random seeds, unit tests, fold-level audits, and numerical checkpoints. The complete pipeline can be run with `python run_all.py`; archived outputs are checked with `python verify_results.py`.

## Results

### Cohort audit established two non-poolable but methodologically complementary contrasts

The LCPM public supplement contained 589 aligned profiles across 676 species-level quantitative features, with 205 CTL, 337 ADE, and 47 CRC participants. The primary CRC-versus-CTL contrast included 252 participants and 336 clean candidate species before prevalence filtering. The MetaCardis public tables contained 1,882 published membership rows but only 1,087 unique participant identifiers. Resolving repeated memberships yielded 994 complete profiles across 729 MGS features. The disjoint primary MetaCardis contrast comprised 303 IHD and 369 MMC participants, with 416 clean bacterial candidates (Figure 1).

The quantitative scales were not interchangeable. LCPM QMP represented estimated species-level cells per gram, whereas MetaCardis QMP represented load-corrected FPKM-derived abundance. Accordingly, all statistical and predictive analyses were conducted within cohort, and no raw abundance values were pooled.

### QMP did not provide a consistent disease-discrimination advantage

In LCPM, mean repeated-cross-validation ROC AUC was 0.659 (SD 0.039) for QMP, 0.653 (SD 0.036) for row-closed abundance, 0.642 (SD 0.036) for minimum-positive CLR, and 0.663 (SD 0.036) for multiplicative-replacement CLR (Figure 2A). The mean paired QMP-minus-row-closed difference was +0.006; its 10 observed repeat-level differences ranged from −0.053 to +0.065. Minimum-positive CLR was 0.017 below QMP on average, whereas multiplicative CLR was 0.004 above QMP, demonstrating sensitivity to zero replacement. At the fixed 0.5 threshold, QMP sensitivity was low (mean 0.204) despite high specificity (0.907), consistent with the 47-versus-205 class imbalance and indicating that these models are not clinically deployable.

In MetaCardis, mean ROC AUC was 0.639 (SD 0.019) for QMP, 0.647 (SD 0.014) for row-closed abundance, 0.643 (SD 0.012) for minimum-positive CLR, and 0.652 (SD 0.013) for multiplicative-replacement CLR (Figure 2B). The mean QMP-minus-row-closed difference reversed direction to −0.008, with observed repeat-level differences from −0.032 to +0.021. QMP-minus-minimum-positive CLR averaged −0.003, and QMP-minus-multiplicative CLR averaged −0.013. Thus, QMP was numerically above row closure in LCPM but below it in MetaCardis. These correlated resampling summaries show no consistent QMP advantage; they do not establish equivalence or non-inferiority.

### QMP and row-closed effects were concordant, but false-discovery-rate calls were representation sensitive

Under the primary pooled, outcome-blind LCPM filter, 93 species were tested in each representation. QMP and row-closed rank-biserial effects were strongly correlated (Pearson r=0.990), with 94.6% directional agreement (Figure 3A). Only *Alistipes onderdonkii* passed FDR<0.05 under all four representations: QMP rank-biserial effect 0.306 (q=0.0406), row-closed effect 0.327 (q=0.0159), minimum-positive CLR effect 0.338 (q=0.0285), and multiplicative-CLR effect 0.352 (q=0.0156).

The source-aligned diagnostic-group-union rule retained 112 species and produced eight QMP and eight row-closed calls, while each CLR rule still produced one call (Figure 4A). The additional seven QMP/row-closed calls therefore depended on an outcome-informed feature universe and are reported only as sensitivity results, not as primary CRC discoveries. None of the LCPM results was adjusted for participant-level BMI, faecal moisture, or calprotectin. The primary result is consequently one unadjusted association plus evidence that filtering decisions materially alter the discovery set; it is not biomarker validation.

In core-adjusted MetaCardis analyses of 404 outcome-blind-filtered MGS features, QMP and row-closed non-zero effects were also strongly correlated (r=0.963), although directional agreement was lower (79.7%; Figure 3B). Despite this concordance, QMP identified six non-zero abundance features at q<0.05 whereas row-closed abundance identified none. The prevalence component identified 10 features; minimum-positive CLR identified 13 and multiplicative-replacement CLR identified 14, with 11 CLR calls shared between the two replacement rules. Representation and zero handling therefore affected thresholded discovery despite broadly concordant QMP and row-closed effect patterns.

### Medication adjustment dominated MetaCardis feature-level conclusions

Adding four medication categories to the core MetaCardis covariates reduced significant prevalence features from 10 to 0, QMP non-zero features from 6 to 0, row-closed non-zero features from 0 to 1, minimum-positive CLR features from 13 to 3, and multiplicative-CLR features from 14 to 0 (Figure 4B). QMP and row-closed effect estimates remained strongly correlated, so the reduction in discoveries was not a wholesale reversal of effect ordering. The two CLR rules also ceased to agree under the expanded specification. These results demonstrate sensitivity to covariate specification but do not identify whether medication is a confounder, mediator, treatment-indication proxy, or marker of disease severity.

### Clinical variables outperformed microbiome-only models in MetaCardis

The estimator-matched MetaCardis clinical-only gradient-boosting model achieved a mean ROC AUC of 0.894 (SD 0.005), substantially higher than microbiome-only AUCs of 0.64–0.65. The row-closed-plus-clinical model reached 0.896 (SD 0.005), an average increase of 0.0018; its 10 repeat-level differences relative to clinical-only ranged from −0.0027 to +0.0095 (Figure 4C). QMP-plus-clinical and minimum-positive-CLR-plus-clinical models achieved AUCs of 0.884 and 0.886, respectively, averaging 0.0099 and 0.0078 below the matched clinical-only model. A logistic-regression clinical baseline produced a similar mean AUC of 0.893. These descriptive internal comparisons show that measured clinical and medication variables dominated discrimination in this contrast, but they do not establish the presence or absence of incremental clinical utility.

### Total microbial load alone showed little disease separation

In LCPM, median total microbial load was 9.16×10¹⁰ cells/g in CTL (interquartile range, 5.29×10¹⁰–1.21×10¹¹) and 8.61×10¹⁰ cells/g in CRC (5.72×10¹⁰–1.21×10¹¹). The Mann–Whitney P value was 0.725, the rank-biserial effect was −0.033, and load alone yielded an AUC of 0.483 (stratified-bootstrap 95% confidence interval, 0.396–0.572; Supplementary Figure 3).

In the MetaCardis load-available subset, median load was 1.10×10¹¹ cells/g in 366 MMC participants (7.38×10¹⁰–1.52×10¹¹) and 1.18×10¹¹ cells/g in 298 IHD participants (8.98×10¹⁰–1.61×10¹¹). The Mann–Whitney P value was 0.0152 and the rank-biserial effect was 0.109, but direct-load discrimination remained weak (AUC 0.555; 0.510–0.599; Supplementary Figure 3). These exploratory findings are consistent with the limited predictive differences between QMP and row-closed abundance, although they do not prove that load is irrelevant to multivariable community structure.

### Representation sensitivity did not reproduce across shared species

Strict harmonization identified 49 exact, one-to-one species shared between the two outcome-blind primary feature sets. The descriptive cross-cohort Pearson and Spearman correlations of standardized QMP-minus-row-closed representation gaps were 0.069 and 0.123, respectively. For QMP-minus-minimum-positive-CLR gaps, the corresponding correlations were 0.180 and 0.209 (Supplementary Figure 1). These weak descriptive correlations do not support consistent species-level representation sensitivity across the two disease settings. Because CRC and IHD are biologically different outcomes and taxa are interdependent, no correlation P values were used and no biological replication claim was made.

## Discussion

This paired reanalysis yielded three main findings. First, QMP did not consistently outperform its row-closed or CLR transformations for internally cross-validated disease discrimination; the QMP-minus-row-closed difference changed direction across cohorts, and multiplicative CLR was numerically highest in both. Second, feature-level conclusions were highly conditional on analysis specification. In LCPM, replacing the source-aligned group-union feature set with pooled outcome-blind filtering reduced QMP and row-closed discoveries from eight to one. In MetaCardis, representation, CLR zero replacement, and the covariate set all altered FDR calls. Third, clinical and medication variables dominated MetaCardis discrimination and association sensitivity. The supported claim is therefore not that one representation is superior, or that predictive performances are equivalent, but that conclusions must be evaluated under explicitly paired preprocessing and sensitivity analyses.

These findings align with previous evidence that microbial load can expose ecological variation obscured by relative profiles [2] and that experimental quantification can reduce compositional and sampling-depth biases under controlled conditions [3]. They also accord with benchmarks showing substantial disagreement among differential-abundance methods and strong effects of confounding [4,5]. Our contribution is deliberately narrower: two public load-informed cohorts were used to quantify how representation, filtering, zero handling, and covariate choice affect association and prediction when participants, feature identities, and cross-validation splits are held constant. This is methodological replication across disease contexts, not a new biological analysis of CRC or IHD.

QMP and row-closed effect estimates remained highly correlated, but this concordance must be interpreted cautiously. Row closure was calculated deterministically from QMP, so the two matrices share the same upstream measurements and technical error. High correlation is partly expected and cannot stand in for a comparison between QMP and an independently measured native sequencing relative profile. The post hoc load results—AUC 0.483 in LCPM and 0.555 in MetaCardis—are consistent with limited load-only discrimination in these contrasts, but they do not establish that microbial load is biologically unimportant or that closure is lossless.

The deterministic construction nevertheless gives this comparison a useful, deliberately narrow interpretation. By holding upstream taxonomic estimates and their measurement error fixed, it isolates the consequences of removing sample-level scale through row-wise closure from differences introduced by independently generated assays. This controlled transformation benchmark asks whether retaining microbial-load scale changes inference or discrimination within these cohorts. It cannot predict how QMP would compare with a native sequencing relative profile, whose library preparation, sampling variation, and taxon-specific measurement errors may differ; that question requires paired independent measurements.

The LCPM filtering sensitivity materially changed the inferential story. Independent filtering can improve multiple-testing efficiency when the filter is unrelated to the test statistic under the null [14], but defining the feature universe from diagnostic-group-specific prevalence uses outcome information. The pooled filter removed that dependence and left one unadjusted species association under all four representations. The eight-call source-aligned result remains useful for reproducing the published workflow but should not be presented as the primary discovery set. More broadly, this result shows that apparent representation effects can be entangled with the feature universe chosen before testing.

CLR also answered a different question from the QMP and non-zero row-closed models because each feature is expressed relative to the geometric mean of retained features. Its dependence on the retained feature set is compounded by the need to replace zeros. Minimum-positive and multiplicative replacement yielded similar LCPM association counts but different MetaCardis call counts and different prediction means. Neither rule can be declared correct from these data. The sensitivity supports reporting the zero-handling rule explicitly and avoiding biological interpretation that depends on only one replacement scheme [1,15].

The MetaCardis medication analysis further shows that quantification does not eliminate clinical confounding. Antihypertensive and lipid-lowering medication use differed sharply between IHD and MMC in the public metadata, and clinical variables alone strongly discriminated the groups. Most feature-level calls disappeared after adding medication indicators, but this does not identify direct drug effects or a deconfounded disease effect. Medication may precede or follow disease, reflect treatment indication, or proxy severity and healthcare context. A causal estimand, prespecified directed acyclic graph, and preferably longitudinal data would be required before interpreting adjustment causally.

The prediction results also require restraint. Microbiome-only AUCs near 0.64–0.66 indicate modest internal discrimination, not a clinically useful diagnostic test. In LCPM, sensitivity at a fixed threshold of 0.5 was approximately 0.20. The estimator-matched MetaCardis comparison showed only a +0.0018 mean AUC difference for row-closed-plus-clinical versus clinical-only, with repeat-level differences in both directions. Because repeated cross-validation reuses participants, these values are descriptive and cannot support equivalence, non-inferiority, or a definitive incremental-utility claim. No external same-disease cohort, calibration analysis, decision curve, or clinically justified threshold was available.

Strengths include one-to-one pairing of representations, explicit auditing of overlapping MetaCardis labels, pooled outcome-blind association filters, fixed classifier specifications, identical participant splits, training-fold-only preprocessing, two CLR zero-replacement rules, FDR control, an estimator-matched clinical comparison, and a reproducible public pipeline with fixed checksums, a patch-level lockfile, unit tests, fold audits, and numerical checkpoints. Keeping cohorts separate prevented inappropriate pooling of cells-per-gram and load-corrected FPKM-derived values. The strict 49-species mapping supplied a transparent descriptive cross-cohort sensitivity analysis without implying biological replication.

Several limitations remain. First, the row-closed matrix is a deterministic transformation of QMP rather than native sequencing relative abundance; shared measurement error limits the strongest comparison. Second, participant-level BMI, faecal moisture, and calprotectin were unavailable for LCPM, leaving its association analysis unadjusted. Third, the primary CRC group was small (n=47), limiting predictive precision and threshold performance. Fourth, the two CLR rules do not exhaust defensible zero-handling methods. Fifth, repeated-cross-validation metrics are correlated descriptive summaries, and no independent test set was available. Sixth, medication covariates lack an identified causal role. Seventh, LCPM and MetaCardis concern different diseases and quantitative scales, precluding a biological meta-analysis. Finally, there was no same-disease external validation, prospective cohort, experimental verification, or independent native-RMP/QMP comparison.

The most informative next study would analyze independently generated native sequencing relative abundance and absolute abundance from flow cytometry, qPCR, or spike-ins, such as the paired technical designs described by Galazzo et al. [13], together with complete covariates and a held-out same-disease cohort. Such a study should prespecify the causal adjustment set, maintain outcome-blind preprocessing, compare justified zero-handling rules, and evaluate calibration and clinically relevant decisions. Those additions would test whether the present methodological sensitivities generalize when quantitative and relative measurements are not deterministically linked.

In conclusion, QMP supplied microbial-load information but did not consistently improve internal disease discrimination in LCPM and MetaCardis. Feature-level conclusions changed with outcome-blind filtering, representation, CLR zero replacement, and covariate specification. Studies with load information should therefore name the estimand represented by each matrix, compare representations on identical participants and splits, report preprocessing sensitivities, and avoid treating QMP, row closure, or CLR as universally superior.

## Conclusion

The current evidence supports a methodological, not causal or biomarker, conclusion: QMP had no consistent predictive advantage, while feature-level results depended materially on filtering, representation, zero handling, and covariate choice. QMP, QMP-derived row closure, and CLR should be compared as distinct estimands under paired, leakage-controlled designs. The results do not establish equivalence and should not be generalized beyond these two cohorts without independent native-profile and same-disease validation.

## Data availability

All input data are publicly available as supplementary workbooks to the original LCPM and MetaCardis publications [6,7]. The analysis pipeline downloads the two source workbooks from Springer Nature and verifies fixed SHA-256 checksums before processing. No controlled-access or participant-level raw sequencing data were redistributed in this project.

## Code availability

All preprocessing, association, prediction, synthesis, verification, and figure-generation code is available at: https://github.com/Joe0908/quantitative-microbiome-compositionality. The repository includes both supported dependency ranges and the exact patch-level reference environment in `requirements-lock.txt`. Before journal submission, the accepted code and manuscript version should be tagged and archived with a permanent DOI (for example, through Zenodo).

## Ethics statement

This study reanalyzed anonymized public supplementary data and involved no new participant recruitment or specimen collection. Ethical approvals and informed-consent procedures for the original cohorts are reported in the source publications [6,7]. Whether the author's institution requires a formal secondary-analysis exemption statement should be confirmed before submission.

## Author contributions

**Template—replace initials after the final author list is confirmed.** Conceptualization: [XX]. Data curation: [XX]. Formal analysis: [XX]. Investigation: [XX]. Methodology: [XX]. Software: [XX]. Validation: [XX]. Visualization: [XX]. Writing—original draft: [XX]. Writing—review and editing: [XX]. Supervision: [XX]. All authors approved the final manuscript and are accountable for the integrity of the work.

## Competing interests

The authors declare no competing interests. **This statement must be confirmed by every author before submission.**

## Funding

[FUNDING SOURCE OR “This research received no specific grant from any funding agency” TO CONFIRM.]

## Acknowledgements

We thank the investigators and participants of the LCPM and MetaCardis studies for making processed data publicly available. [ADDITIONAL ACKNOWLEDGEMENTS TO CONFIRM.]

## Generative-AI disclosure

OpenAI ChatGPT/Codex was used to assist with code organization, figure preparation, manuscript structuring, and language editing under author direction. All numerical claims were linked to auditable analysis outputs, and the human authors remain responsible for verification, interpretation, authorship, and the submitted text. This statement should be adapted to the target journal's current policy.

## Figure legends

### Figure 1. Study design and paired abundance representations

Public LCPM and MetaCardis quantitative microbiome matrices were audited separately and transformed into matched published quantitative (QMP), QMP-derived row-closed, minimum-positive CLR, and multiplicative-replacement CLR representations. LCPM provided species-level cells-per-gram values; MetaCardis provided a microbial-load-corrected FPKM-derived index. Primary contrasts were CRC versus CTL in LCPM (n=252) and IHD versus MMC in MetaCardis (n=672). Pooled outcome-blind filtering retained 93 LCPM species and 404 MetaCardis MGS features for association analysis. Inference, prediction, and robustness analyses were performed within cohort; raw values and disease effects were not pooled. Row-closed abundance is a deterministic QMP transformation, not a native sequencing relative profile. CRC, colorectal cancer; CTL, lesion-free control; IHD, ischaemic heart disease; MMC, metabolically matched control; MGS, metagenomic species; CLR, centered log-ratio.

### Figure 2. Repeated-cross-validation discrimination across abundance representations

Repeat-level out-of-fold ROC AUCs for QMP, QMP-derived row-closed, minimum-positive CLR, and multiplicative-replacement CLR representations in (A) LCPM CRC versus CTL and (B) MetaCardis IHD versus MMC. Each point represents one of 10 repeats of stratified five-fold cross-validation; gray lines connect results from the same repeat. Large points and error bars show the mean and sample standard deviation across repeats. All representations used identical participant splits and fold-local pooled feature filtering; minimum-positive CLR references were fitted on training folds. The dashed horizontal line indicates AUC=0.5. Repeat-level summaries are descriptive because participants and training sets recur across repeats; no equivalence or non-inferiority inference is implied.

### Figure 3. Concordance of QMP and row-closed feature effects

(A) LCPM rank-biserial effects for 93 pooled-outcome-blind-filtered species in CRC versus CTL. (B) Core-adjusted MetaCardis non-zero log-abundance disease coefficients for 404 pooled-outcome-blind-filtered MGS features in IHD versus MMC. Blue points were significant under QMP or row-closed abundance at Benjamini–Hochberg q<0.05; gray points were not. Dashed lines indicate equality. Correlations describe effect-pattern concordance; row closure is deterministic from QMP and the measurement units are not independent.

### Figure 4. Filtering, covariate, and clinical-model sensitivity

(A) LCPM FDR-significant call counts under the pooled outcome-blind 93-species filter and the source-aligned, diagnostic-group-union 112-species filter. (B) MetaCardis call counts under core and core-plus-medication covariate specifications for prevalence, QMP non-zero, row-closed non-zero, minimum-positive CLR, and multiplicative-replacement CLR components. Medication-expanded results are an alternative-covariate sensitivity, not causal drug-effect estimates. (C) Mean repeated-cross-validation ROC AUC for estimator-matched MetaCardis clinical-only and microbiome-plus-clinical histogram-gradient-boosting models; error bars show sample standard deviations across 10 correlated repeats. Panel C is descriptive and does not establish incremental clinical utility.

### Supplementary Figure 1. Cross-cohort representation sensitivity across exact shared species

Within-cohort standardized representation gaps for 49 strict one-to-one species shared by the outcome-blind primary LCPM and MetaCardis feature sets. (A) QMP-minus-row-closed gaps. (B) QMP-minus-minimum-positive-CLR gaps. Cross-cohort correlations are descriptive because taxa are interdependent. They do not compare or pool CRC and IHD disease effects.

### Supplementary Figure 2. LCPM prevalence-threshold and test sensitivity

(A) Number of CTL–ADE–CRC features with Benjamini–Hochberg q<0.05 at 1%, 5%, and 10% prevalence thresholds under the source-aligned diagnostic-group-union rule, using QMP and row-closed Kruskal–Wallis tests and minimum-positive CLR one-way analysis of variance or Kruskal–Wallis tests. (B) Number of features retained at each prevalence threshold. This outcome-informed analysis is supplementary and does not define the primary CRC feature set.

### Supplementary Figure 3. Total microbial-load sensitivity

ROC AUC of total microbial load alone, with 95% stratified percentile-bootstrap confidence intervals from 10,000 resamples. The dashed line indicates AUC=0.5. LCPM included 205 CTL and 47 CRC participants (n=252); MetaCardis included 366 MMC and 298 IHD participants with available load (n=664). This analysis was post hoc and was performed separately by cohort.

## Supplementary Methods

### Source-file verification

The LCPM and MetaCardis supplementary workbooks were downloaded by HTTPS and accepted only if their SHA-256 hashes matched `2d9c0fa807fbcd85f97beee292d24551920d33bc76435c4ea578f5d90cc10282` and `ce68279db4ce3c0c29a244ef6fb5ff712dcee3e6c3b5b33250174368b0c74248`, respectively. A checksum mismatch stopped the pipeline.

### Fold-level leakage audit

For each of 50 folds per cohort, the pipeline stored train/test sample counts, class counts, selected-feature count, and minimum and maximum minimum-positive CLR replacement values. No outcome label was used by the pooled prevalence filter; no held-out abundance or clinical value was used to fit feature filtering, BMI imputation, or minimum-positive replacement vectors. Multiplicative replacement was row-wise after the fold-local feature set was fixed.

### Repeated-cross-validation interpretation

Five held-out folds were concatenated within each repeat so every participant contributed one out-of-fold probability. Performance was calculated on this full repeat-level vector. The 10 values were summarized by their mean, sample standard deviation, and observed minimum-to-maximum range. No t interval, Wilcoxon test, P value, or q value was calculated across repeats because the estimates share participants and overlapping training data. The analysis was not designed as an equivalence or non-inferiority study.

### Exact-species harmonization

LCPM feature labels were converted from dot-separated to space-separated binomials. MetaCardis features were eligible only when their species label was a strict two-token binomial and did not contain generic “sp.” or “bacterium” labels. Species occurring more than once in either primary outcome-blind source set were excluded. The resulting 49 mappings were inspected as one-to-one links before standardized representation-gap analysis.

## Supplementary tables and files plan

1. **Supplementary Table S1:** Cohort and source-file audit, inclusion/exclusion flow, public group overlaps, and measurement units.
2. **Supplementary Table S2:** Primary-group descriptive statistics, missingness, and medication prevalence for MetaCardis; diagnosis counts and microbial load for LCPM.
3. **Supplementary Table S3:** Complete LCPM CRC-versus-CTL QMP, row-closed, minimum-positive CLR, and multiplicative-CLR results under the 93-species primary and 112-species source-aligned filters.
4. **Supplementary Table S4:** Complete MetaCardis prevalence, QMP non-zero, row-closed non-zero, and both CLR results for all 404 features under all three covariate specifications.
5. **Supplementary Table S5:** Repeat- and fold-level cross-validation metrics, descriptive paired differences, and leakage-audit fields.
6. **Supplementary Table S6:** Exact mappings and standardized effects for the 49 shared species.
7. **Supplementary Table S7:** Total microbial-load sensitivity statistics and bootstrap intervals.
8. **Supplementary Figure S1:** Shared-species representation-gap comparison.
9. **Supplementary Figure S2:** LCPM prevalence-threshold and CLR-test sensitivity.
10. **Supplementary Figure S3:** Total microbial-load sensitivity.
11. **Supplementary reporting file:** Completed STORMS checklist with manuscript page references.

## References

1. Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. Microbiome datasets are compositional: and this is not optional. *Front Microbiol.* 2017;8:2224. https://doi.org/10.3389/fmicb.2017.02224
2. Vandeputte D, Kathagen G, D'Hoe K, et al. Quantitative microbiome profiling links gut community variation to microbial load. *Nature.* 2017;551:507–511. https://doi.org/10.1038/nature24460
3. Lloréns-Rico V, Vieira-Silva S, Gonçalves PJ, et al. Benchmarking microbiome transformations favors experimental quantitative approaches to address compositionality and sampling depth biases. *Nat Commun.* 2021;12:3562. https://doi.org/10.1038/s41467-021-23821-6
4. Nearing JT, Douglas GM, Hayes MG, et al. Microbiome differential abundance methods produce different results across 38 datasets. *Nat Commun.* 2022;13:342. https://doi.org/10.1038/s41467-022-28034-z
5. Wirbel J, Essex M, Forslund SK, et al. A realistic benchmark for differential abundance testing and confounder adjustment in human microbiome studies. *Genome Biol.* 2024;25:247. https://doi.org/10.1186/s13059-024-03390-9
6. Tito RY, et al. Microbiome confounders and quantitative profiling challenge predicted microbial targets in colorectal cancer development. *Nat Med.* 2024;30:1339–1348. https://doi.org/10.1038/s41591-024-02963-2
7. Fromentin S, Forslund SK, Chechi K, et al. Microbiome and metabolome features of the cardiometabolic disease spectrum. *Nat Med.* 2022;28:303–314. https://doi.org/10.1038/s41591-022-01688-4
8. Topçuoğlu BD, Lesniak NA, Ruffin MT IV, Wiens J, Schloss PD. A framework for effective application of machine learning to microbiome-based classification problems. *mBio.* 2020;11:e00434-20. https://doi.org/10.1128/mBio.00434-20
9. Garach Vélez I, Ortuño Guzmán FM, Rojas Ruiz I, Herrera Maldonado LJ. Exploring the role of normalization and feature selection in microbiome disease classification pipelines. *GigaScience.* 2025;14:giaf096. https://doi.org/10.1093/gigascience/giaf096. Funding correction: *GigaScience.* 2026;15:giag050. https://doi.org/10.1093/gigascience/giag050
10. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B.* 1995;57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
11. Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python. *J Mach Learn Res.* 2011;12:2825–2830. https://jmlr.org/papers/v12/pedregosa11a.html
12. Mirzayi C, Renson A, Genomic Standards Consortium, et al. Reporting guidelines for human microbiome research: the STORMS checklist. *Nat Med.* 2021;27:1885–1892. https://doi.org/10.1038/s41591-021-01552-x
13. Galazzo G, Tedjo DI, Wintjens DSJ, et al. How to count our microbes? The effect of different quantitative microbiome profiling approaches. *Front Cell Infect Microbiol.* 2020;10:403. https://doi.org/10.3389/fcimb.2020.00403
14. Bourgon R, Gentleman R, Huber W. Independent filtering increases detection power for high-throughput experiments. *Proc Natl Acad Sci USA.* 2010;107:9546–9551. https://doi.org/10.1073/pnas.0914005107
15. Martín-Fernández JA, Barceló-Vidal C, Pawlowsky-Glahn V. Dealing with zeros and missing values in compositional data sets using nonparametric imputation. *Math Geol.* 2003;35:253–278. https://doi.org/10.1023/A:1023866030544
