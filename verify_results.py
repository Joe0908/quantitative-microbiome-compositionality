"""Verify a completed run against the documented reference checkpoints."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from compositionality.config import DEFAULT_OUTPUT_DIR, PROJECT_ROOT


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _call_map(frame: pd.DataFrame) -> dict[str, int]:
    return dict(
        frame[["component", "significant_q_lt_0_05"]]
        .itertuples(index=False, name=None)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--auc-tolerance", type=float, default=5e-4)
    args = parser.parse_args()

    expected = json.loads((PROJECT_ROOT / "expected_results.json").read_text())
    cv = pd.read_csv(args.output_dir / "synthesis" / "microbiome_cv_summary.csv")
    all_cv = pd.read_csv(args.output_dir / "prediction" / "cv_summary.csv")
    for cohort, models in expected["microbiome_roc_auc_means"].items():
        for model, target in models.items():
            observed = float(
                cv.loc[
                    (cv.cohort == cohort) & (cv.model == model), "roc_auc_mean"
                ].iloc[0]
            )
            require(
                abs(observed - target) <= args.auc_tolerance,
                f"{cohort} {model} AUC: expected {target}, observed {observed}",
            )

    for cohort, target in expected["clr_multiplicative_roc_auc_means"].items():
        observed = float(
            all_cv.loc[
                (all_cv.cohort == cohort)
                & (all_cv.model == "CLR (multiplicative)"),
                "roc_auc_mean",
            ].iloc[0]
        )
        require(
            abs(observed - target) <= args.auc_tolerance,
            f"{cohort} multiplicative CLR AUC: expected {target}, observed {observed}",
        )

    for model, target in expected["optional_metacardis_roc_auc_means"].items():
        row = all_cv.loc[(all_cv.cohort == "MetaCardis") & (all_cv.model == model)]
        if row.empty:  # valid after ``--microbiome-only``
            continue
        observed = float(row["roc_auc_mean"].iloc[0])
        require(
            abs(observed - target) <= args.auc_tolerance,
            f"MetaCardis {model} AUC: expected {target}, observed {observed}",
        )

    da = pd.read_csv(
        args.output_dir / "synthesis" / "differential_association_counts.csv"
    )
    for specification, expected_key in [
        ("pooled_outcome_blind", "LCPM_pooled_outcome_blind"),
        ("source_aligned_group_union", "LCPM_source_aligned_group_union"),
    ]:
        current = da.loc[
            (da.cohort == "LCPM")
            & (da.filter_specification == specification)
        ]
        require(
            _call_map(current)
            == expected["differential_calls_q_lt_0_05"][expected_key],
            f"{expected_key} differential calls differ",
        )
        require(
            set(current["tested_features"])
            == {expected["tested_features"][expected_key]},
            f"{expected_key} tested-feature count differs",
        )

    lcpm_results = pd.read_csv(
        args.output_dir / "associations" / "lcpm_crc_vs_ctl_associations.csv"
    )
    pooled_significant = lcpm_results.loc[
        (lcpm_results["filter_specification"] == "pooled_outcome_blind")
        & lcpm_results["significant_q_lt_0_05"],
        "feature",
    ]
    require(
        sorted(pooled_significant.unique())
        == expected["LCPM_pooled_primary_significant_species"],
        "LCPM pooled primary significant-species set differs",
    )

    for adjustment, expected_key in [
        ("core", "MetaCardis_core"),
        ("core_plus_medications", "MetaCardis_core_plus_medications"),
    ]:
        current = da.loc[
            (da.cohort == "MetaCardis") & (da.adjustment == adjustment)
        ]
        require(
            _call_map(current)
            == expected["differential_calls_q_lt_0_05"][expected_key],
            f"{expected_key} differential calls differ",
        )
        require(
            set(current["tested_features"])
            == {expected["tested_features"]["MetaCardis_pooled_outcome_blind"]},
            f"{expected_key} tested-feature count differs",
        )

    descriptive = pd.read_csv(
        args.output_dir / "prediction" / "cv_paired_descriptive_differences.csv"
    )
    forbidden = {
        "ci95_low",
        "ci95_high",
        "paired_wilcoxon_p",
        "paired_wilcoxon_q_bh",
    }
    require(
        forbidden.isdisjoint(descriptive.columns),
        "Repeated-CV output still contains invalid inferential columns",
    )
    require(
        not any("ci95" in column for column in all_cv.columns),
        "Repeated-CV summary still contains inferential confidence intervals",
    )

    shared = pd.read_csv(args.output_dir / "synthesis" / "shared_exact_species.csv")
    require(len(shared) == expected["shared_exact_species"], "Shared-species count differs")

    headline = json.loads(
        (args.output_dir / "synthesis" / "headline_results.json").read_text()
    )
    observed_correlations = headline["qmp_row_closed_effect_correlations"]
    require(
        abs(
            observed_correlations["LCPM"]
            - expected["qmp_row_closed_effect_correlations"]["LCPM"]
        )
        < 1e-9,
        "LCPM QMP/row-closed effect correlation differs",
    )
    require(
        abs(
            observed_correlations["MetaCardis"]
            - expected["qmp_row_closed_effect_correlations"]["MetaCardis_core"]
        )
        < 1e-9,
        "MetaCardis QMP/row-closed effect correlation differs",
    )
    print("All reference checkpoints passed.")


if __name__ == "__main__":
    main()
