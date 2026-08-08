"""Reproducible post hoc sensitivity analysis of total microbial load.

This analysis is intentionally separate from the prespecified QMP/row-closed/CLR
pipeline. It asks whether total microbial load alone separates the two primary
disease contrasts and provides a mechanistic check on the observed similarity
of representation-specific classifiers.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score


SEED = 531
BOOTSTRAP_REPLICATES = 10_000
ROOT = Path(__file__).resolve().parents[2]
PIPELINE = ROOT
OUTPUT_DIR = Path(__file__).resolve().parent


def stratified_auc_interval(
    negative: np.ndarray,
    positive: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float]:
    aucs = np.empty(BOOTSTRAP_REPLICATES, dtype=float)
    y = np.r_[np.zeros(len(negative), dtype=int), np.ones(len(positive), dtype=int)]
    for index in range(BOOTSTRAP_REPLICATES):
        sampled_negative = rng.choice(negative, len(negative), replace=True)
        sampled_positive = rng.choice(positive, len(positive), replace=True)
        aucs[index] = roc_auc_score(y, np.r_[sampled_negative, sampled_positive])
    return tuple(np.quantile(aucs, [0.025, 0.975]))


def summarize(
    cohort: str,
    contrast: str,
    negative_label: str,
    positive_label: str,
    negative: np.ndarray,
    positive: np.ndarray,
    rng: np.random.Generator,
) -> dict[str, object]:
    test = stats.mannwhitneyu(negative, positive, alternative="two-sided", method="asymptotic")
    # scipy reports U for the first input; positive-minus-negative rank-biserial
    # therefore reverses the first-input orientation.
    rank_biserial = 1.0 - 2.0 * float(test.statistic) / (len(negative) * len(positive))
    y = np.r_[np.zeros(len(negative), dtype=int), np.ones(len(positive), dtype=int)]
    values = np.r_[negative, positive]
    auc = float(roc_auc_score(y, values))
    low, high = stratified_auc_interval(negative, positive, rng)
    return {
        "cohort": cohort,
        "contrast": contrast,
        "negative_group": negative_label,
        "positive_group": positive_label,
        "n_negative": int(len(negative)),
        "n_positive": int(len(positive)),
        "negative_median_cells_per_g": float(np.median(negative)),
        "negative_q1_cells_per_g": float(np.quantile(negative, 0.25)),
        "negative_q3_cells_per_g": float(np.quantile(negative, 0.75)),
        "positive_median_cells_per_g": float(np.median(positive)),
        "positive_q1_cells_per_g": float(np.quantile(positive, 0.25)),
        "positive_q3_cells_per_g": float(np.quantile(positive, 0.75)),
        "mann_whitney_p": float(test.pvalue),
        "rank_biserial_positive_minus_negative": rank_biserial,
        "direct_load_auc": auc,
        "bootstrap_auc_ci95_low": float(low),
        "bootstrap_auc_ci95_high": float(high),
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "seed": SEED,
        "analysis_status": "post hoc exploratory sensitivity analysis",
    }


def main() -> None:
    rng = np.random.default_rng(SEED)

    lcpm = pd.read_csv(PIPELINE / "data/processed/lcpm/metadata.csv.gz")
    ctl = lcpm.loc[lcpm["diagnosis"].eq("CTL"), "total_load_cells_per_g"].to_numpy(float)
    crc = lcpm.loc[lcpm["diagnosis"].eq("CRC"), "total_load_cells_per_g"].to_numpy(float)

    meta = pd.read_csv(PIPELINE / "data/processed/metacardis/metadata.csv.gz")
    meta = meta.loc[
        meta["quantitative_profile_available"].astype(bool)
        & meta["microbial_load_available"].astype(bool)
        & (meta["ihd_member"].astype(bool) | meta["mmc_member"].astype(bool))
    ]
    mmc = meta.loc[meta["mmc_member"].astype(bool), "microbial_load_cells_per_g"].to_numpy(float)
    ihd = meta.loc[meta["ihd_member"].astype(bool), "microbial_load_cells_per_g"].to_numpy(float)

    rows = [
        summarize("LCPM", "CRC vs CTL", "CTL", "CRC", ctl, crc, rng),
        summarize("MetaCardis", "IHD372 vs MMC372", "MMC372", "IHD372", mmc, ihd, rng),
    ]
    frame = pd.DataFrame(rows)
    frame.to_csv(OUTPUT_DIR / "load_sensitivity_summary.csv", index=False)
    (OUTPUT_DIR / "load_sensitivity_summary.json").write_text(
        json.dumps(rows, indent=2), encoding="utf-8"
    )
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
