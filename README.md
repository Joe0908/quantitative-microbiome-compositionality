# Quantitative versus compositional microbiome pipeline

This repository is the complete, separated-code version of the Galazzo/LCPM +
MetaCardis project. It replaces the old relative-abundance-only dataset with two
published quantitative resources, derives paired RMP and CLR views from the
same participants, and tests how representation changes prediction and
taxon-level inference.

The raw FASTQ files are **not required** for this project. Reprocessing reads
would introduce a much larger bioinformatics project and would not, by itself,
recover absolute abundance. The public quantitative supplementary matrices are
the appropriate starting data.

## What the two datasets contribute

| Dataset | Primary comparison | Profiles | Quantitative value | Role |
|---|---:|---:|---|---|
| Galazzo/LCPM | CRC vs CTL | 47 vs 205 | Species cells per gram | True quantitative benchmark |
| MetaCardis | IHD372 vs MMC372 | 303 vs 369 | Microbial-load-corrected MGS abundance index | Independent methodological replication |

MetaCardis QMP is **not cells per gram**. Therefore, diseases, raw QMP values,
and disease-effect estimates are never pooled across cohorts.

## Quick start

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python run_all.py
```

The first run downloads and SHA-256-verifies both public supplementary
workbooks. To use workbooks already placed in `data/raw/`, run:

```bash
python run_all.py --skip-download
```

To reproduce only the primary QMP/RMP/CLR comparison and skip the optional
MetaCardis clinical prediction baselines:

```bash
python run_all.py --microbiome-only
```

## Code parts

Each analytical part is an independent module and also runs from
`run_all.py`.

| Part | Module | Content |
|---:|---|---|
| 01 | `part01_download_data.py` | Download, cache, and checksum both public workbooks |
| 02 | `part02_prepare_lcpm.py` | Build LCPM metadata, QMP cells/g, derived RMP, and taxonomy |
| 03 | `part03_prepare_metacardis.py` | Stream the large workbook, audit overlapping labels, and build QMP/RMP |
| 04 | `part04_lcpm_associations.py` | CRC-vs-CTL QMP/RMP/CLR tests plus prevalence sensitivity |
| 05 | `part05_metacardis_associations.py` | Prevalence, non-zero QMP/RMP, and CLR regression models |
| 06 | `part06_train_models.py` | Leakage-controlled repeated cross-validation |
| 07 | `part07_synthesize.py` | Within-cohort comparison and guarded exact-species synthesis |

Examples for running one part:

```bash
python -m compositionality.part02_prepare_lcpm
python -m compositionality.part05_metacardis_associations
python -m compositionality.part06_train_models --microbiome-only
python -m compositionality.part07_synthesize
```

## Expected primary checkpoints

| Cohort | QMP AUC | RMP AUC | CLR AUC | Differential calls |
|---|---:|---:|---:|---|
| Galazzo/LCPM | 0.658983 | 0.652870 | 0.641692 | 8 / 8 / 1 |
| MetaCardis, core adjusted | 0.639364 | 0.647203 | 0.642858 | 6 / 0 / 13 |

For MetaCardis, core-adjusted prevalence has 10 calls. Adding the four broad
medication categories gives 0 prevalence, 0 QMP, 1 RMP, and 3 CLR calls.
There are 51 one-to-one exact species shared across the final tested sets.

Small last-decimal differences can occur across supported versions of
scikit-learn or statsmodels. `expected_results.json` records the reference
environment checkpoints.

## Main outputs

The repository does not redistribute the large public source workbooks or
generated result tables. Running the pipeline creates the following local
directories; `expected_results.json` contains the verified checkpoints:

```text
data/
  raw/
  processed/lcpm/
  processed/metacardis/
outputs/
  associations/
  prediction/
  synthesis/
```

Important result files are:

- `outputs/prediction/cv_summary.csv`
- `outputs/prediction/cv_paired_comparisons.csv`
- `outputs/associations/lcpm_crc_vs_ctl_qmp_rmp_clr.csv`
- `outputs/associations/metacardis_ihd_vs_mmc_hurdle_qmp_rmp_clr.csv`
- `outputs/synthesis/headline_results.json`
- `outputs/synthesis/shared_exact_species.csv`

## Reproducibility and interpretation rules

- Repeated stratified CV is 5 folds × 10 repeats with seed 531.
- The same split and fold-specific feature set is used for QMP, RMP, and CLR.
- Prevalence filtering and CLR zero replacement are fitted on training data
  only.
- No hyperparameter tuning, calibration, external validation, or threshold
  selection is performed.
- Repeat-level confidence intervals and paired Wilcoxon tests are exploratory,
  because repeated-CV estimates are correlated.
- Quantitative abundance does not remove confounding. The medication analysis
  is a sensitivity check, not proof that medication is a pure confounder.
- Shared-species results describe representation sensitivity. They do not imply
  that CRC and IHD have the same microbiome biology.

See `PROJECT_OUTLINE.md` for the complete data and training specification.

## Public sources

- Galazzo/LCPM study: <https://www.nature.com/articles/s41591-024-02963-2>
- MetaCardis study: <https://www.nature.com/articles/s41591-022-01688-4>

The exact supplementary URLs and SHA-256 values are fixed in
`compositionality/config.py`.
