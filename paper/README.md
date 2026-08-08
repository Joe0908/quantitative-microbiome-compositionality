# Manuscript and submission package

This directory contains the reviewer-driven manuscript revision matched to the verified analysis outputs on this branch.

## Main files

- `manuscript.md` — editable manuscript source.
- `research_audit_and_submission.md` — research audit, reviewer critique, revision matrix, journal positioning, cover letter, and submission checklist.
- `analysis/` — scripts and audited outputs used for the post hoc microbial-load analysis and all manuscript figures.
- `build_documents.py` — reproducible DOCX builder.

Rendered figures and DOCX files are generated artifacts and are intentionally not versioned. Running the commands below recreates `paper/figures/`, `paper/Quantitative_Microbiome_Manuscript_Draft.docx`, and `paper/Research_Audit_and_Submission_Package.docx`.

## Regenerate figures

Run the main pipeline first so that `outputs/` and `data/processed/` exist, then run:

```bash
python paper/analysis/compute_load_sensitivity.py
python paper/analysis/make_figures.py
```

## Regenerate documents

Install `python-docx`, then run:

```bash
python paper/build_documents.py
```

## Scientific status

The revision fixes the statistical and reporting issues identified during external review:

- repeated-cross-validation results are descriptive only; no repeat-level P values, q values, or confidence intervals are reported;
- association filtering is pooled and outcome blind in the primary analysis;
- row-closed abundance is explicitly identified as a deterministic QMP transformation, not native sequencing relative abundance;
- minimum-positive and multiplicative CLR zero replacement are both evaluated;
- the MetaCardis clinical-only comparison uses an estimator-matched histogram-gradient-boosting model;
- medication-expanded models are treated as non-causal sensitivity analyses.

The manuscript is not submission-ready until author metadata, funding/conflict statements, the formatted supplement, a completed STORMS checklist, and a permanent release DOI are supplied.
