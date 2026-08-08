# Research Audit and Submission Package

## Executive decision

**Current status: C — Publishable after moderate revision for a realistic methods/replication journal.**

The study can support a defensible paper if it is positioned as a paired, cross-disease methodological reanalysis. It cannot support claims of a new CRC biomarker panel, IHD mechanism, causal microbial effect, universal superiority of QMP, or external biological replication. An ambitious microbial-ecology journal would probably require a third cohort with independently measured native relative and quantitative profiles, or same-disease external validation.

## 1. Research Audit

### A. Research question

> When the same stool microbiome profiles are represented as published quantitative abundance, QMP-derived row-closed composition, or centered log-ratio abundance, how much do disease discrimination and feature-level conclusions change?

**Scientific value:** High. The question addresses compositionality, microbial load, differential association, and machine-learning reproducibility.

**Novelty:** Moderate. Prior work has established the compositionality problem and benchmarked transformations. The incremental contribution is a rigorously paired real-data benchmark across two load-informed public cohorts with fold-local preprocessing and a transparent confounding sensitivity analysis.

**Answerable with current data:** Yes, for internal methodological comparison. No, for causal inference, clinical utility, or same-disease external validation.

**Publication suitability:** Suitable for a technically focused, soundness-oriented microbiology or bioinformatics journal after the remaining reproducibility and reporting items are completed.

### B. Central claim

> Across LCPM and MetaCardis, QMP showed no consistent internal discrimination advantage, while feature-level conclusions depended materially on outcome-blind filtering, abundance representation, CLR zero replacement, and covariate specification.

This claim is fully supported by the archived results. It deliberately avoids claiming that one representation is correct, that the diseases share microbial effects, or that the models are clinically useful.

### Title candidates

1. **Selected:** A paired benchmark of quantitative, row-closed, and log-ratio gut microbiome profiles in two public cohorts
2. Outcome-blind filtering and abundance representation shape microbiome inference in two load-informed cohorts
3. Quantitative microbiome profiling shows no consistent discrimination advantage in paired analyses of LCPM and MetaCardis
4. Filtering, zero handling, and abundance representation alter feature discovery in two quantitative microbiome cohorts
5. Paired quantitative and compositional representations yield unstable inferential but similar predictive rankings across two cohorts

### C. Novelty audit

#### What is already known

1. Sequencing-derived relative microbiome profiles are compositional and can produce misleading marginal interpretations.
2. Experimental quantitative profiling can reveal microbial-load-related ecological variation.
3. Differential-abundance methods often disagree on the same data.
4. Confounding can dominate microbiome disease associations.
5. Normalization, feature selection, and classifier choice interact in disease-classification pipelines.
6. The source LCPM publication already compared quantitative, relative, and CLR analyses and examined CRC confounders.

#### What remains insufficiently characterized

Real-data evidence is still limited on whether the same abundance-representation choice affects feature-level inference and multivariable discrimination to the same degree when participant splits, feature identities, and leakage controls are explicitly paired.

#### Gap filled by this study

The project benchmarks prediction and feature-level inference in two public load-informed cohorts, audits MetaCardis participant overlap, uses pooled outcome-blind association filters and training-fold-only prediction preprocessing, tests two CLR zero-replacement rules, and separates method replication from biological replication.

#### Increment beyond prior work

- Paired QMP, row-closed, and two CLR prediction specifications under identical repeated-cross-validation splits.
- Pooled outcome-blind association filtering plus fold-local prediction filtering and CLR preprocessing.
- Explicit separation of MetaCardis prevalence and non-zero abundance components.
- Core-versus-medication sensitivity analysis.
- Estimator-matched MetaCardis clinical-only versus clinical-plus-microbiome comparison.
- Strict cross-cohort mapping of 49 species without pooling disease effects.
- Public end-to-end code with fixed input checksums, a patch-level lockfile, tests, and numerical checkpoints.

#### Novelty category

| Category | Assessment |
|---|---|
| Conceptual | Low–moderate; the central distinction between inference and prediction is useful but not wholly new. |
| Methodological | Moderate; the paired validation design and audit are the strongest contribution. |
| Empirical | Moderate; two real load-informed cohorts support the pattern. |
| Dataset-related | Low; both datasets were already published. |
| Translational | Low; no clinical validation or deployment analysis was performed. |

### D. Data audit

#### Data structure

| Domain | LCPM | MetaCardis |
|---|---|---|
| Public source | Supplementary tables to Tito et al. 2024 | Supplementary tables to Fromentin et al. 2022 |
| Full audited participants | 589 | 1,087 unique IDs from 1,882 membership rows |
| Complete quantitative profiles | 589 | 994 |
| Primary contrast | 205 CTL vs 47 CRC | 369 MMC vs 303 IHD |
| Quantitative scale | Species-level cells/g | Load-corrected FPKM-derived index |
| Source features | 676 species | 729 MGS features |
| Clean candidates | 336 | 416 |
| Primary association features | 93 pooled outcome-blind; 112 source-aligned sensitivity | 404 pooled outcome-blind |
| Key covariates | Not public at participant level | Age, BMI, sex, nationality, diabetes, four drug groups |
| External same-disease validation | No | No |

#### Critical problems

1. **Derived relative matrix, not native sequencing RMP.** Row-closed abundance is a deterministic transform of QMP. It shares measurement error with QMP and cannot test disagreement between independently measured sequencing RMP and absolute abundance. This is fatal only if the paper claims an experimental QMP-versus-RMP comparison. The manuscript now uses “QMP-derived row-closed abundance” throughout.
2. **Novelty overlap with the source LCPM paper.** The original publication already evaluated QMP, relative, and CLR associations and confounding. A biomarker-focused CRC paper would be redundant. The current manuscript instead centers the paired inferential-versus-predictive benchmark.
3. **No same-disease external validation.** MetaCardis concerns IHD, not CRC. It validates the methodological question only. It cannot validate LCPM species or CRC prediction.
4. **LCPM confounder data are unavailable.** Participant-level BMI, stool moisture, and calprotectin values could not be adjusted. Significant LCPM species cannot be promoted as deconfounded CRC biomarkers.

#### Important problems

1. The CRC class is small (n=47) and the fixed-threshold sensitivity is approximately 0.20. Precision and clinical usefulness are limited.
2. MetaCardis IHD and MMC differ substantially in medications, sex, and country. The high clinical AUC is partly a study-design and treatment signature, not necessarily disease biology.
3. Medication adjustment may control confounding, block mediators, or encode treatment indication. A causal estimand and adjustment set were not prespecified.
4. CLR remains dependent on zero handling. Minimum-positive and multiplicative replacement were tested, but neither is a gold standard.
5. Repeated-cross-validation values are correlated. The revision removed t intervals, Wilcoxon tests, P values, and q values across repeats; means, SDs, and observed ranges are descriptive only.
6. MetaCardis MGS features are not necessarily unique species. Discovery counts must be called “features,” not unique taxa.
7. Batch effects could not be re-estimated from the public processed tables. Source-study processing may have addressed them, but residual batch cannot be ruled out.

#### Minor problems

1. A permanent software/data release DOI is not yet available.
2. A completed STORMS checklist and final formatted supplementary tables are not yet assembled.
3. Author names, affiliations, funding, contribution initials, conflict statements, and the institutional secondary-analysis position remain to be confirmed.

#### Missing analyses and explicit status

- **External same-disease validation:** “This result has not yet been calculated.” A suitable independent CRC cohort with quantitative abundance is required.
- **Calibration curves and calibration slope/intercept:** “This result has not yet been calculated.” Participant-level out-of-fold probabilities must be retained or regenerated.
- **Decision-curve analysis and a clinically justified threshold:** “This result has not yet been calculated.” A target clinical use and threshold must first be defined.
- **Native sequencing RMP versus independently measured QMP:** “This result has not yet been calculated.” A dataset such as Galazzo 2020 or another paired technical cohort is needed.
- **Formal prospective power calculation:** “This result has not yet been calculated.” The current secondary analysis was not designed prospectively; power should be based on a defined primary estimand for a future validation study.

## 2. Reconstructed scientific narrative

### Logic

**Background** → Relative microbiome data are compositional; QMP restores load information.

**Knowledge gap** → It is unclear whether representation changes feature-level inference and multivariable discrimination to the same extent in real load-informed cohorts.

**Research question** → Compare published QMP, QMP-derived row closure, and CLR on identical participants and features.

**Working expectation** → QMP will not be universally superior; inferential and predictive sensitivity will be quantified without assuming one must exceed the other.

**Study design** → LCPM CRC-vs-CTL and MetaCardis IHD-vs-MMC, analyzed separately, with paired association and leakage-controlled prediction workflows.

**Main result 1** → QMP, row-closed, and CLR AUCs were similar within cohort; the QMP advantage reversed direction.

**Why needed** → Directly tests the predictive component of the hypothesis.

**Main result 2** → QMP and row-closed effects were highly correlated, but LCPM calls depended strongly on the filter and MetaCardis calls depended on representation and CLR replacement.

**Why needed** → Shows that effect-pattern stability and discovery stability are not equivalent.

**Main result 3** → Medication adjustment removed most MetaCardis calls, while clinical variables dominated discrimination.

**Why needed** → Tests whether quantification eliminates confounding; it does not.

**Robustness** → Total-load AUCs were weak, filtering and CLR replacement changed discovery counts, and representation gaps across 49 shared species correlated weakly across cohorts.

**Why needed** → Defines when QMP may have limited added predictive information and prevents false cross-disease replication claims.

**Interpretation** → Representations encode different estimands. QMP is informative but not universally predictive; row closure and CLR remain necessary sensitivity analyses.

**Limitations** → Derived row closure, missing LCPM covariates, no same-disease external validation, zero-replacement assumptions, correlated CV summaries.

**Conclusion** → No representation was universally superior; filtering, zero handling, and covariates materially shape feature-level conclusions, while the internal predictive ranking also varied by cohort and CLR rule.

### Reviewer-feedback triage completed in this revision

| Advice | Decision | Evidence/action |
|---|---|---|
| Remove inferential tests across repeated CV runs | Accepted | Removed all repeat-level t intervals, Wilcoxon tests, P values, and q values; retained mean, SD, and observed range only. |
| Describe row closure as native relative abundance | Rejected | It is deterministic from QMP and is now named “QMP-derived row-closed abundance” throughout. |
| Add outcome-blind LCPM filtering | Accepted with correction | Used pooled prevalence across all CRC/CTL participants. CTL-only filtering was not used because it still conditions on outcome labels. |
| Add CLR zero-replacement sensitivity | Accepted | Added multiplicative replacement with δ=1/p² alongside minimum-positive replacement. |
| Match the clinical and combined estimators | Accepted | Added a clinical-only histogram-gradient-boosting model on identical folds while retaining logistic regression as context. |
| Treat medication adjustment as a negative control | Rejected | Medication has no established causal role here; the expanded model is described only as an alternative-covariate sensitivity. |
| Add a causal DAG immediately | Deferred | A DAG without a prespecified estimand and adequate temporal variables would create false causal precision. Causal ambiguity is stated explicitly. |
| Package as `pip install` for mSystems | Not required | The limiting issue for mSystems is novelty/external validation, not install syntax. A one-command reproducible workflow and lockfile are sufficient for this revision. |

## 3. Figure-driven manuscript plan

| Figure | Scientific purpose | Input/sample size | Statistical display | Axes/groups | Main takeaway | Results subsection |
|---|---|---|---|---|---|---|
| Figure 1 | Show cohorts, measurement units, participant flow, representations, and analysis boundaries | LCPM full n=589/primary n=252; MetaCardis 1,087 IDs/994 profiles/primary n=672 | Descriptive workflow | Cohort → representation → inference/prediction/robustness | The analyses are paired within cohort and not pooled across disease or scale | Cohort audit |
| Figure 2 | Test whether QMP consistently improves discrimination and whether CLR zero handling changes ranking | LCPM n=252; MetaCardis n=672; 10 repeated-CV AUCs per model | Paired repeat AUC trajectories; mean±SD | QMP, row-closed, minimum-positive CLR, multiplicative CLR | No consistent QMP advantage; zero handling changes CLR means | Disease discrimination |
| Figure 3 | Separate effect concordance from discovery concordance | 93 LCPM species; 404 MetaCardis features | QMP effect vs row-closed effect; Pearson r; FDR colors | x=QMP effect; y=row-closed effect | Effects correlate strongly even when FDR calls differ | Feature-level effects |
| Figure 4A | Quantify filtering sensitivity | LCPM 93 pooled vs 112 source-aligned species | FDR call counts | Representation by filter | Outcome-informed group-union filtering expands QMP/row-closed calls | LCPM filtering sensitivity |
| Figure 4B | Show covariate and zero-handling sensitivity | MetaCardis 404 features | FDR call counts | Component by covariate set | Medication-expanded and CLR-replacement results differ sharply | Medication sensitivity |
| Figure 4C | Place microbiome discrimination in estimator-matched clinical context | MetaCardis n=672 | Repeated-CV AUC mean±SD | Clinical HGB and three combined HGB models | Clinical variables dominate this study contrast | Clinical context |
| Supplementary Figure 1 | Describe cross-cohort representation sensitivity | 49 exact shared species | Gap-vs-gap scatter; Pearson/Spearman | LCPM gap vs MetaCardis gap | Species-level gaps correlate weakly across diseases | Cross-cohort synthesis |
| Supplementary Figure 2 | Test threshold and CLR-statistic sensitivity | LCPM CTL/ADE/CRC; 1%, 5%, 10% | Call counts and feature-set size | Prevalence threshold | Discovery depends on filter and CLR test | Robustness |
| Supplementary Figure 3 | Test whether total load alone separates disease groups | LCPM n=252; MetaCardis load subset n=664 | AUC with 10,000-replicate stratified-bootstrap CI | Cohort | Load is non-discriminatory or weakly discriminatory | Microbial-load sensitivity |

## 4. Statistical review

| Analysis | Question | Input | Method | Main assumptions | Multiple testing | Effect/CI | Required robustness or alternative |
|---|---|---|---|---|---|---|---|
| LCPM CRC vs CTL | Do feature distributions differ? | 205 CTL, 47 CRC; 93 pooled primary features | Two-sided Mann–Whitney U | Independent participants; ordinal/rank-comparable observations | BH within filter/representation | Rank-biserial effect; q value | Covariate adjustment if participant data become available |
| LCPM CLR | Does feature log-ratio position differ? | Row-closed 93-feature primary matrix | Minimum-positive or multiplicative replacement + CLR + Mann–Whitney | Retained feature set and replacement define the reference | BH by replacement family | Rank-biserial effect | Both planned rules completed; no gold-standard replacement claim |
| LCPM three-group sensitivity | Are CTL/ADE/CRC distributions different? | Full n=589 | Kruskal–Wallis; CLR ANOVA and KW | Independent groups; ANOVA residual assumptions for CLR | BH per threshold/representation/test | Test statistic and q | Keep both CLR tests; report that author supplement labels conflict with F-distribution behavior only in reproducibility notes |
| MetaCardis prevalence | Does probability of detection differ? | 672 participants, 404 features | Binomial GLM | Correct link; independent participants; no severe separation | BH per adjustment/component | Log-odds coefficient | Penalized/Firth sensitivity only if rare-feature claims become central |
| MetaCardis non-zero abundance | Among carriers, does abundance differ? | Feature-specific positive subsets | log abundance OLS + HC3 t inference | Linear mean structure; positive-subset selection is understood; robust SE adequate | BH per adjustment/component | Log-scale coefficient and SE | A unified two-part model or MaAsLin 3 sensitivity could strengthen results |
| MetaCardis CLR | Does relative log-ratio position differ? | All 672 participants | CLR OLS + HC3 | Zero replacement and geometric-mean reference are appropriate | BH | CLR coefficient and SE | Alternative zero replacement; country-stratified sensitivity |
| QMP/row-closed effect stability | Are effect patterns similar? | Matched feature effects | Pearson r; sign agreement; Jaccard calls | Feature estimates treated descriptively; correlations not independent biological observations | Not primary testing | r and agreement | Spearman r and bootstrap over features could be supplementary |
| Prediction | Does representation change discrimination? | LCPM n=252; MetaCardis n=672 | Fixed HGB, repeated stratified 5×10 CV | Participants are independent; repeated metrics are correlated; no leakage | No inferential testing across repeats | AUC mean±SD and observed range | External validation is essential; no equivalence claim |
| Clinical context | Does microbiome add to clinical variables? | MetaCardis n=672 | Estimator-matched HGB clinical-only and combined models; logistic context | Same folds/estimator for matched arm; repeat metrics correlated | No inferential testing across repeats | Descriptive AUC difference/range | External validation, calibration, and a clinical target would be required for utility claims |
| Total microbial load | Does load alone separate groups? | LCPM n=252; Meta n=664 | MWU, rank-biserial, direct-score AUC, stratified bootstrap | Load is measured comparably within cohort; direct monotonic score is appropriate | Two exploratory tests; no confirmatory family | AUC 95% bootstrap CI | Adjusted load model; nonlinear load terms; treat as post hoc |
| Shared species | Is representation sensitivity species-level similar? | 49 one-to-one species | Within-cohort z scores; Pearson/Spearman gap correlations | Mapping is exact; disease effects are not pooled; taxa are interdependent | Descriptive; no feature-correlation P value | Correlation only | Third cohort; do not meta-analyze CRC and IHD effects |

## 5. Reviewer Mode

### Reviewer 1 — Domain expert

#### Major comments

1. The biological novelty is limited because both cohorts and quantitative matrices were previously published, and the LCPM source paper already analyzed QMP, relative abundance, CLR, and confounding.
2. The paper must not present the eight LCPM species as newly discovered CRC biomarkers. Without BMI, moisture, and calprotectin adjustment, the analysis is less biologically controlled than the source publication.
3. MetaCardis is a different disease setting and cannot validate CRC biology. The manuscript should consistently say “methodological replication.”
4. Row-closed values are derived from QMP, so their high correlation is partly built in. A native relative profile would be more biologically informative.
5. The explanation for why feature discovery changes more than AUC should be framed as an interpretation, not a demonstrated mechanism.

#### Minor comments

1. Define MMC at first mention and explain why it is the comparator.
2. Use “MGS features,” not “species,” for MetaCardis discovery counts.
3. Clarify that total-load analysis was post hoc.
4. Avoid “validation” in the title.
5. Provide a complete source-data provenance table.

### Reviewer 2 — Statistical reviewer

#### Major comments

1. Repeated-CV folds and repeats are correlated, so the t intervals and Wilcoxon P values should not be presented as conventional inferential confidence intervals.
2. LCPM association filtering uses outcome groups and should be repeated with an outcome-blind prevalence filter or explicitly marked exploratory.
3. The clinical-only logistic model and combined gradient-boosting models are not estimator matched; incremental utility is not identified.
4. CLR zero replacement is arbitrary. At least one alternative replacement sensitivity should be considered if CLR differences are central.
5. MetaCardis medication adjustment lacks a causal model. Drug covariates may be confounders, mediators, or treatment-indication proxies.
6. The small CRC group limits performance precision and threshold metrics. No external validation is available.
7. Patch-level environment versions were not retained.

#### Minor comments

1. Provide fold-level selected-feature counts and replacement ranges.
2. Report class prevalence alongside average precision.
3. Explain why no hyperparameter tuning was performed.
4. State the FDR family for every q value.
5. Add status counts for non-estimable MetaCardis fits.

### Reviewer 3 — Highly skeptical reviewer

#### Strongest argument for rejection

> The main comparison is partly deterministic because the “relative” matrix is calculated from QMP, the central biological datasets have already been extensively analyzed by their original authors, and the second cohort concerns a different disease; therefore the manuscript may be viewed as a technically careful but incrementally novel reanalysis without genuine external validation.

#### Major comments

1. The study does not establish that quantitative profiling improves biological truth, because no gold standard or independent native relative measurement is analyzed.
2. Similar AUCs may reflect an underpowered or weak classifier rather than representation equivalence.
3. Different feature-level calls may be driven by sparsity, zero replacement, or test estimands rather than compositionality itself.
4. Clinical confounding in MetaCardis is so strong that the IHD contrast may primarily classify treatment and country.
5. The study lacks a prospective validation cohort and formal replication of any disease-specific association.
6. The current repository nomenclature incorrectly implies that Galazzo was analyzed.

#### Minor comments

1. The title should not imply all quantitative microbiome datasets.
2. Provide a permanent release DOI.
3. Reduce mechanistic discussion of individual taxa.
4. Add an explicit table of what can and cannot be inferred.

## 6. Revision Mode

| Reviewer concern | Severity | Can current data address it? | Required analysis/action | Manuscript change |
|---|---|---|---|---|
| Row-closed matrix is not native RMP | Critical | Partly | No new analysis can remove deterministic linkage; obtain paired native RMP/QMP for a stronger paper | Use “QMP-derived row-closed”; state shared-error limitation everywhere |
| LCPM novelty overlaps source paper | Critical | Yes through framing | Center inference-vs-prediction benchmark, not CRC biomarkers | Title, abstract, introduction, discussion already reframed |
| No same-disease external validation | Critical for ambitious/clinical claim | No | Add independent quantitative CRC or IHD cohort | Explicitly limit to methodological replication |
| Missing LCPM confounders | Critical for biomarker claim | No with public supplement | Request data or remove biomarker claim | Treat LCPM feature findings as unadjusted reanalysis |
| Outcome-informed LCPM association filter | Important | Resolved with current data | Pooled outcome-blind primary filter; source group-union retained as sensitivity | Primary calls changed from 8/8/1/1 to 1/1/1/1 |
| Correlated repeated-CV inference | Important | Resolved in reporting | Removed all t intervals, Wilcoxon tests, P values, and q values across repeats | Mean, SD, and actual range only; no equivalence claim |
| Mixed clinical estimators | Important | Resolved with current data | Added clinical-only HGB under the same folds/settings | Logistic baseline retained as context only |
| CLR zero replacement | Important | Partly resolved | Added multiplicative replacement sensitivity | Report both rules; no gold-standard claim |
| Medication causal ambiguity | Important | Partly | Do not draw a causal DAG without a defined estimand; state ambiguity | Alternative-covariate sensitivity only, not causal control |
| Exact environment missing | Important | Resolved | Added and tested `requirements-lock.txt` | Report exact patch versions |
| “Galazzo/LCPM” code label | Important for integrity | Resolved | Replaced public-facing label and regenerated outputs | Use “LCPM” consistently |
| Supplement and STORMS checklist missing | Important for submission | Yes | Assemble tables S1–S7 and completed checklist | Add page-linked supplement |
| Authorship/ethics/funding placeholders | Blocking administrative | Yes, author input needed | Confirm details and institutional policy | Replace every bracketed field |
| No archive DOI | Minor–important | Yes | Tag release and archive in Zenodo | Update code availability and citation |

### Prioritized revision checklist

1. Confirm author list, affiliations, corresponding author, funding, conflicts, contributions, and institutional secondary-analysis requirements.
2. Assemble Supplementary Tables S1–S7 and complete the STORMS checklist.
3. Tag the reviewed GitHub release and create a permanent DOI.
4. Recheck journal-specific word limits, figure rules, and generative-AI policy at submission.
5. For mSystems or ISME Communications, add a truly independent native-RMP/QMP dataset or same-disease external validation.
6. Do not add more classifiers, thresholds, or exploratory taxa without a claim they directly answer.

## 7. Internal consistency audit

### Sample numbers

| Location/analysis | Required number |
|---|---:|
| LCPM full public profiles | 589 |
| LCPM CTL/ADE/CRC | 205/337/47 |
| LCPM primary CRC-vs-CTL | 252 |
| MetaCardis published rows / unique IDs | 1,882 / 1,087 |
| MetaCardis complete profiles | 994 |
| MetaCardis primary MMC/IHD | 369/303; total 672 |
| MetaCardis load subset MMC/IHD | 366/298; total 664 |
| LCPM primary/source-aligned association features | 93/112 |
| MetaCardis primary association features | 404 |
| Shared exact species | 49 |

### Primary statistics

| Result | Required value |
|---|---|
| LCPM AUC QMP/row-closed/CLR-min/CLR-mult | 0.659/0.653/0.642/0.663 |
| LCPM QMP-minus-row-closed | +0.006; observed repeat range −0.053 to +0.065; no P/q |
| MetaCardis AUC QMP/row-closed/CLR-min/CLR-mult | 0.639/0.647/0.643/0.652 |
| MetaCardis QMP-minus-row-closed | −0.008; observed repeat range −0.032 to +0.021; no P/q |
| LCPM primary calls QMP/row/CLR-min/CLR-mult | 1/1/1/1 |
| LCPM source-aligned sensitivity calls | 8/8/1/1 |
| MetaCardis core calls prevalence/QMP/row/CLR-min/CLR-mult | 10/6/0/13/14 |
| MetaCardis medication calls prevalence/QMP/row/CLR-min/CLR-mult | 0/0/1/3/0 |
| QMP–row-closed effect correlations | 0.990 LCPM; 0.963 MetaCardis |
| Clinical HGB and row-closed+clinical HGB AUC | 0.894 and 0.896 |
| Load-only AUC LCPM/MetaCardis | 0.483/0.555 |

### Terminology rules

- Use **LCPM**, never “Galazzo/LCPM.”
- Use **row-closed abundance** or **QMP-derived row-closed composition**, not native RMP.
- Use **MGS feature** for MetaCardis feature counts.
- Use **methodological replication**, not external biological validation.
- Use **association**, **suggests**, or **is consistent with**, not causation.
- Use **observed repeat range**, not confidence interval, P value, q value, equivalence, or non-inferiority language for repeated-CV summaries.

## 8. Publication-readiness assessment

| Dimension | Score (/10) | Rationale |
|---|---:|---|
| Scientific question | 8.0 | Clear, relevant, and appropriately narrow after reframing |
| Novelty | 5.5 | Moderate paired-design contribution, but substantial prior and source-paper overlap |
| Dataset | 7.0 | Two strong public cohorts; different diseases/scales and no same-disease external cohort |
| Methods | 8.0 | Outcome-blind filtering, paired folds, two CLR rules, and matched clinical estimator; deterministic row closure remains limiting |
| Statistics | 8.0 | Repeat-level pseudo-inference removed; FDR families explicit; causal-adjustment and external-validation limitations remain |
| Results | 8.0 | Revised conclusions track verified outputs, including the substantive 8-to-1 LCPM filtering change |
| Biological interpretation | 7.5 | Carefully bounded; intentionally not a mechanistic paper |
| Figures | 8.0 | Four focused main figures and three targeted supplementary sensitivities |
| Writing | 8.0 | Clear and conservative; requires target-journal formatting and author edits |
| Reproducibility | 9.0 | Public code, checksums, seed, lockfile, tests, fold audits, and numerical verification; DOI still pending |
| Overall publishability | 7.2 | Realistic after moderate revision; not ready for an ambitious journal without independent validation |

**Current status: C — Publishable after moderate revision.**

## 9. Journal positioning

### Ambitious

#### 1. mSystems

**Scope fit:** Official scope includes computational microbiology and synthesis across large multidimensional datasets: https://journals.asm.org/journal/msystems/scope

**Why it fits:** The paired representation benchmark, confounding analysis, and public workflow address microbial-systems methodology.

**Likely weakness:** Current work is a reanalysis of two published cohorts and lacks independent native relative abundance or same-disease validation.

**Needed before submission:** Add a paired native-RMP/QMP dataset or external same-disease validation; the locked environment is complete, but a permanent archive remains necessary.

#### 2. ISME Communications

**Scope fit:** The journal welcomes methodological and computational studies that represent major advances in microbial ecology: https://academic.oup.com/ismecommun/pages/about

**Why it fits:** The study addresses how microbial-community measurement representation changes inference.

**Likely weakness:** Editors may judge the advance too narrow because row closure is deterministic and the diseases differ.

**Needed before submission:** Demonstrate generalization with an independent quantitative method and strengthen ecological interpretation.

### Realistic

#### 1. Microbiology Spectrum

**Scope fit:** Broad basic, applied, and clinical microbial sciences, with decisions focused on technical quality and community usefulness; its FAQ explicitly welcomes replication and negative results: https://journals.asm.org/journal/spectrum/scope and https://journals.asm.org/journal/spectrum/faq

**Why it fits:** This is a technically careful replication/robustness study whose important result is the absence of a universal QMP predictive advantage.

**Likely weakness:** The editor may still request a clearer independent contribution beyond source analyses.

**Needed before submission:** Code naming and the environment lock are complete. Finish author metadata, supplement, STORMS checklist, and permanent release DOI.

#### 2. GigaScience

**Scope fit:** Its official criteria emphasize open, FAIR data, reproducibility, usability, and utility for life-science and biomedical research: https://academic.oup.com/gigascience/pages/About

**Why it fits:** The checksummed public inputs, end-to-end workflow, lockfile, audits, and negative/robustness result align closely with open-science priorities.

**Likely weakness:** The study reuses published processed tables and does not introduce a new algorithm or large new data resource.

**Needed before submission:** Deposit the versioned workflow and reusable outputs with persistent identifiers, complete the journal's reporting checklist, and make the benchmark utility clearer than the biological novelty.

#### 3. Scientific Reports

**Scope fit:** Broad soundness-based research with public-data reanalyses and computational studies.

**Why it fits:** The analysis is quantitative, reproducible, and appropriately non-causal.

**Likely weakness:** Reviewers may request more biological interpretation or independent validation.

**Needed before submission:** Final supplement, exact environment, and careful title/abstract framing.

### Safe

#### 1. PLOS ONE

**Scope fit:** Technical validity and ethical reporting rather than predicted impact.

**Why it fits:** The paper provides a reproducible negative/robustness result from public data.

**Likely weakness:** Methods and software provenance must be exceptionally complete.

**Needed before submission:** Permanent archive DOI, environment lock, supplement, and completed reporting checklist.

#### 2. Frontiers in Microbiology — relevant microbiome or bioinformatics section

**Scope fit:** Broad microbiome and computational microbiology research.

**Why it fits:** Compositionality and quantitative profiling are directly relevant.

**Likely weakness:** Section fit and article-processing costs should be checked; reviewers may ask for another quantitative cohort.

**Needed before submission:** Same core reporting revisions; new analyses are optional if claims remain modest.

## 10. Cover letter draft — Microbiology Spectrum

Dear Editor,

Please consider our manuscript, “A paired benchmark of quantitative, row-closed, and log-ratio gut microbiome profiles in two public cohorts,” for publication as a Research Article in *Microbiology Spectrum*.

Microbiome sequencing data are compositional, whereas quantitative profiling retains information about microbial load. We performed a paired reanalysis of two public load-informed gut-microbiome datasets: LCPM colorectal cancer versus lesion-free controls and MetaCardis ischaemic heart disease versus metabolically matched controls. Published quantitative abundance, its deterministically derived row-closed composition, and centered log-ratio abundance under two zero-replacement rules were compared using pooled outcome-blind association filters, identical prediction splits, training-fold-only preprocessing, and fixed classifiers.

The central result is deliberately nuanced. QMP did not show a consistent predictive advantage: the descriptive QMP-minus-row-closed mean AUC difference was +0.006 in LCPM and −0.008 in MetaCardis. More importantly, pooled outcome-blind filtering reduced the primary LCPM QMP/row-closed discoveries from eight under the source-aligned group-union rule to one. CLR findings varied with zero replacement, medication-expanded adjustment removed most MetaCardis discoveries, and clinical variables dominated discrimination. These results show that quantitative measurement is informative but does not eliminate sensitivity to preprocessing, estimand choice, or clinical context.

The study is suited to *Microbiology Spectrum* because its current scope explicitly considers useful replication studies, negative results, reanalyses of large datasets, and methodological work. The manuscript does not claim causal effects, clinical utility, equivalence, or biological replication across diseases. All input tables are public, and the complete pipeline, exact environment lock, checksums, seeds, fold-level audits, and numerical checkpoints are available at https://github.com/Joe0908/quantitative-microbiome-compositionality. A permanent release DOI will be added before submission.

This manuscript is not under consideration elsewhere. All authors have approved the submission and declare [NO COMPETING INTERESTS / DETAILS TO CONFIRM]. No new human participants or specimens were involved; the analysis used anonymized public supplementary data from studies with ethics procedures reported in their original publications.

Thank you for your consideration.

Sincerely,

[CORRESPONDING AUTHOR]

[AFFILIATION]

[EMAIL]

## 11. Suggested reviewers

These suggestions are based on directly relevant methodological expertise and publicly verifiable institutional profiles. They must be conflict-checked against all authors, supervisors, recent collaborators, and shared institutions before use.

1. **Jakob Wirbel, PhD — Helmholtz Centre for Infection Research (HZI).** Expertise: realistic differential-abundance benchmarking, confounder adjustment, microbiome machine learning. Official profile: https://www.helmholtz-hzi.de/en/persons/dr-jakob-wirbel/
2. **Jacob T. Nearing, PhD — Harvard T.H. Chan School of Public Health / Broad Institute (current affiliation should be reconfirmed at submission).** Expertise: differential-abundance method disagreement and microbiome bioinformatics. Profile: https://nearinj.github.io/
3. **Ignacio Garach Vélez — University of Granada.** Expertise: normalization and feature selection in microbiome disease-classification pipelines. Official profile: https://www.ugr.es/en/staff/ignacio-garach-velez
4. **Chloe Mirzayi — CUNY Graduate School of Public Health and Health Policy (current title should be reconfirmed).** Expertise: microbiome reporting, causal inference, and reproducibility; STORMS lead author. Profile: https://cunyisph.org/team/chloe-mirzayi/

Do not suggest investigators from the LCPM or MetaCardis author teams because their direct involvement in the source datasets could create perceived conflicts.

## 12. Submission handoff checklist

- [ ] Replace all author, affiliation, email, funding, contribution, conflict, and acknowledgement placeholders.
- [ ] Obtain approval from every author and confirm authorship order.
- [ ] Confirm whether institutional ethics/exemption wording is required.
- [x] Rename public-facing “Galazzo/LCPM” labels to “LCPM” and regenerate outputs.
- [x] Create a patch-level environment lockfile and rerun the revised analyses.
- [x] Run unit tests and `verify_results.py` successfully.
- [ ] Assemble Supplementary Tables S1–S7.
- [ ] Complete the STORMS checklist with manuscript page references.
- [x] Add pooled outcome-blind filtering, multiplicative CLR replacement, and estimator-matched clinical sensitivity analyses.
- [ ] Tag a GitHub release and mint a Zenodo DOI.
- [ ] Recheck every abstract number against the final tables.
- [ ] Recheck every figure label, legend, and sample number.
- [ ] Confirm all references and journal formatting on the submission date.
- [ ] Check reviewer conflicts and current affiliations.
- [ ] Adapt the generative-AI disclosure to the target journal's current policy.
