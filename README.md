# Quantitative versus compositional microbiome pipeline

This repository contains the complete analysis for a methodological comparison
of quantitative microbiome profiling (QMP), QMP-derived row-closed abundance,
and centered log-ratio (CLR) abundance in the public LCPM and MetaCardis
cohorts.

The row-closed matrix is calculated deterministically from each published QMP
matrix. It is **not** an independently measured, native sequencing relative-
abundance matrix. QMP and row-closed values therefore share upstream
measurement error. The project evaluates the consequences of mathematical
representation, not two independent laboratory assays.

Raw FASTQ files are not required. The analysis starts from the quantitative
supplementary matrices released with the two studies and verifies them by
SHA-256 checksum.

## Cohorts and scope

| Dataset | Primary comparison | Profiles | Published quantitative value | Role |
|---|---:|---:|---|---|
| LCPM | CRC vs CTL | 47 vs 205 | Species cells per gram | Quantitative CRC benchmark |
| MetaCardis | IHD372 vs MMC372 | 303 vs 369 | Load-corrected MGS abundance index | Cross-disease methodological replication |

MetaCardis QMP is not cells per gram. Diseases, raw abundance values, and
disease-effect estimates are never pooled across cohorts. The cross-cohort
analysis is methodological replication, not same-disease biological
validation.

## Quick start

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python run_all.py
```

`requirements-lock.txt` records the exact Python 3.12.13 reference
environment used for the numerical checkpoints. Use it instead of
`requirements.txt` when exact environment reconstruction is required.

The first run downloads and checksum-verifies both public supplementary
workbooks. To use workbooks already placed in `data/raw/`, run:

```bash
python run_all.py --skip-download
```

To skip the optional MetaCardis clinical and combined prediction models:

```bash
python run_all.py --microbiome-only
```

Verify a completed run and execute the unit tests:

```bash
python verify_results.py
python -m pytest -q
```

## Code parts

| Part | Module | Content |
|---:|---|---|
| 01 | `part01_download_data.py` | Download, cache, and checksum both workbooks |
| 02 | `part02_prepare_lcpm.py` | Build LCPM metadata, QMP cells/g, row-closed abundance, and taxonomy |
| 03 | `part03_prepare_metacardis.py` | Audit overlapping labels and build MetaCardis QMP/row-closed data |
| 04 | `part04_lcpm_associations.py` | Outcome-blind CRC–CTL associations and source-aligned sensitivity |
| 05 | `part05_metacardis_associations.py` | Prevalence, positive-abundance, and CLR association models |
| 06 | `part06_train_models.py` | Leakage-controlled repeated cross-validation |
| 07 | `part07_synthesize.py` | Guarded within- and cross-cohort synthesis |

Each part is independently runnable, for example:

```bash
python -m compositionality.part04_lcpm_associations
python -m compositionality.part05_metacardis_associations
python -m compositionality.part06_train_models --microbiome-only
python -m compositionality.part07_synthesize
```

## Revised reference checkpoints

### Prediction

| Cohort | QMP | Row-closed | CLR, minimum-positive | CLR, multiplicative |
|---|---:|---:|---:|---:|
| LCPM | 0.658983 | 0.652870 | 0.641692 | 0.662605 |
| MetaCardis | 0.639364 | 0.647203 | 0.642858 | 0.652397 |

These are means across ten repeated out-of-fold prediction vectors. Repeats
reuse participants and are correlated; the pipeline reports mean, sample SD,
and the literal repeat minimum–maximum only. It does not assign confidence
intervals or inferential P values to repeated-CV differences, and the analysis
is not an equivalence or non-inferiority test.

Using the same HistGradientBoosting estimator and folds in MetaCardis, the
clinical-only AUC is 0.893944. The corresponding combined AUCs are 0.884081
(QMP), 0.895759 (row-closed), and 0.886138 (minimum-positive CLR).

### Feature-level inference

- LCPM primary pooled outcome-blind filter: 93 taxa and 1/1/1/1 FDR calls for
  QMP, row-closed, minimum-positive CLR, and multiplicative CLR.
- LCPM source-aligned group-union sensitivity: 112 taxa and 8/8/1/1 calls.
  This analysis uses diagnosis groups to define the tested set and is not the
  primary inferential result.
- MetaCardis pooled outcome-blind filter: 404 taxa. Core-adjusted calls are 10
  prevalence, 6 QMP, 0 row-closed, 13 minimum-positive CLR, and 14
  multiplicative CLR.
- Adding four medication categories as an alternative adjustment set yields
  0, 0, 1, 3, and 0 calls, respectively. This is a confounding sensitivity
  analysis, not a causal estimate of medication effects.

The revised exact-species synthesis contains 49 one-to-one matches.

## Main outputs

Running the pipeline creates `data/raw/`, `data/processed/`, and `outputs/`.
Generated data and result tables are intentionally not committed.

Important result files include:

- `outputs/associations/lcpm_crc_vs_ctl_associations.csv`
- `outputs/associations/metacardis_ihd_vs_mmc_hurdle_qmp_row_closed_clr.csv`
- `outputs/prediction/cv_summary.csv`
- `outputs/prediction/cv_paired_descriptive_differences.csv`
- `outputs/synthesis/lcpm_outcome_blind_filter_sensitivity.csv`
- `outputs/synthesis/clr_zero_replacement_cv_summary.csv`
- `outputs/synthesis/headline_results.json`
- `outputs/synthesis/shared_exact_species.csv`

## Reproducibility and interpretation rules

- Repeated stratified CV uses 5 folds × 10 repeats with seed 531.
- All abundance representations use the same split and fold-specific taxa.
- Prediction-time prevalence filtering, minimum-positive CLR replacement, and
  BMI imputation are fitted inside each training fold.
- Multiplicative CLR uses `delta = 1 / p²` after closing the selected taxa.
- No hyperparameter tuning, calibration, decision-curve analysis, external
  validation, or threshold selection is performed.
- LCPM taxon tests are unadjusted because the required participant-level BMI,
  stool-moisture, and calprotectin covariates are not public.
- Quantitative measurement does not remove clinical confounding.

See `PROJECT_OUTLINE.md` for the full training and statistical specification.

## Public sources

- LCPM study: <https://www.nature.com/articles/s41591-024-02963-2>
- MetaCardis study: <https://www.nature.com/articles/s41591-022-01688-4>

The exact supplementary URLs and SHA-256 values are fixed in
`compositionality/config.py`.
