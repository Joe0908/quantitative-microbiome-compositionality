"""Part 07 — synthesize the paired LCPM and MetaCardis analyses.

Raw abundances and disease effects are never pooled.  Cross-cohort work is
limited to exact one-to-one species labels and standardized, within-cohort
representation gaps.
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .common import (
    BINOMIAL_DOT,
    BINOMIAL_SPACE,
    ensure_dir,
    paired_repeat_descriptions,
    require,
    write_json,
    zscore,
)
from .config import DEFAULT_OUTPUT_DIR


EXPECTED_AUCS = {
    ("LCPM", "QMP"): 0.6589828749351323,
    ("LCPM", "Row-closed"): 0.6528697457187339,
    ("LCPM", "CLR"): 0.6416917488323820,
    ("MetaCardis", "QMP"): 0.6393642616294150,
    ("MetaCardis", "Row-closed"): 0.6472027690573936,
    ("MetaCardis", "CLR"): 0.6428577817131307,
}

PRIMARY_REPRESENTATIONS = ["QMP", "Row-closed", "CLR"]


def _jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def _stability_rows(
    cohort: str,
    contrast: str,
    adjustment: str,
    table: pd.DataFrame,
    feature_column: str,
    representation_column: str,
) -> list[dict[str, object]]:
    effect = table.pivot(
        index=feature_column, columns=representation_column, values="effect"
    )
    calls = {
        representation: set(
            table.loc[
                (table[representation_column] == representation)
                & table["significant_q_lt_0_05"].astype(bool),
                feature_column,
            ]
        )
        for representation in PRIMARY_REPRESENTATIONS
    }
    rows = []
    for left, right in combinations(PRIMARY_REPRESENTATIONS, 2):
        paired = effect[[left, right]].dropna()
        left_calls, right_calls = calls[left], calls[right]
        rows.append(
            {
                "cohort": cohort,
                "contrast": contrast,
                "adjustment": adjustment,
                "representation_a": left,
                "representation_b": right,
                "calls_a": len(left_calls),
                "calls_b": len(right_calls),
                "call_overlap": len(left_calls & right_calls),
                "call_jaccard": _jaccard(left_calls, right_calls),
                "effect_pearson_r": float(paired[left].corr(paired[right])),
                "direction_agreement": float(
                    (np.sign(paired[left]) == np.sign(paired[right])).mean()
                ),
                "estimable_features": len(paired),
            }
        )
    return rows


def _meta_representation(value: str) -> str | None:
    return {
        "qmp_nonzero": "QMP",
        "row_closed_nonzero": "Row-closed",
        "clr_minimum_positive": "CLR",
    }.get(value)


def _shared_species(
    lcpm: pd.DataFrame, meta: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    lcpm_effect = lcpm.pivot(
        index="feature", columns="representation", values="effect"
    )[PRIMARY_REPRESENTATIONS]
    for representation in PRIMARY_REPRESENTATIONS:
        lcpm_effect[f"{representation}_z"] = zscore(lcpm_effect[representation])
    lcpm_effect = lcpm_effect.reset_index()
    lcpm_effect["species"] = lcpm_effect["feature"].where(
        lcpm_effect["feature"].map(
            lambda value: bool(BINOMIAL_DOT.fullmatch(str(value)))
        )
    )
    lcpm_effect["species"] = lcpm_effect["species"].str.replace(
        ".", " ", regex=False
    )
    lcpm_counts = lcpm_effect["species"].value_counts(dropna=True)
    lcpm_one = lcpm_effect.loc[
        lcpm_effect["species"].map(lcpm_counts).eq(1).fillna(False)
    ].copy()

    meta_core = meta.loc[
        (meta["adjustment"] == "core")
        & meta["component"].isin(
            ["qmp_nonzero", "row_closed_nonzero", "clr_minimum_positive"]
        )
    ].copy()
    meta_core["representation"] = meta_core["component"].map(_meta_representation)
    meta_effect = meta_core.pivot(
        index=["matrix_column", "feature_id", "species"],
        columns="representation",
        values="effect",
    ).reset_index()
    for representation in PRIMARY_REPRESENTATIONS:
        meta_effect[f"{representation}_z"] = zscore(meta_effect[representation])
    strict = meta_effect["species"].map(
        lambda value: bool(BINOMIAL_SPACE.fullmatch(str(value)))
    )
    conservative_cross_cohort = ~meta_effect["species"].astype(str).str.contains(
        r"\bsp\.|bacterium", case=False, regex=True
    )
    meta_strict = meta_effect.loc[strict & conservative_cross_cohort].copy()
    meta_counts = meta_strict["species"].value_counts()
    duplicated_labels = meta_counts.index[meta_counts > 1]
    meta_one = meta_strict.loc[~meta_strict["species"].isin(duplicated_labels)].copy()

    shared = lcpm_one.merge(
        meta_one,
        on="species",
        how="inner",
        suffixes=("_lcpm", "_metacardis"),
        validate="one_to_one",
    )
    require(bool(len(shared)), "No one-to-one exact species were shared across cohorts")
    for cohort_suffix in ["lcpm", "metacardis"]:
        shared[f"qmp_minus_row_closed_z_gap_{cohort_suffix}"] = (
            shared[f"QMP_z_{cohort_suffix}"]
            - shared[f"Row-closed_z_{cohort_suffix}"]
        )
        shared[f"qmp_minus_clr_z_gap_{cohort_suffix}"] = (
            shared[f"QMP_z_{cohort_suffix}"] - shared[f"CLR_z_{cohort_suffix}"]
        )

    summary_rows = []
    for name, lcpm_column, meta_column in [
        (
            "QMP minus row-closed standardized gap",
            "qmp_minus_row_closed_z_gap_lcpm",
            "qmp_minus_row_closed_z_gap_metacardis",
        ),
        (
            "QMP minus CLR standardized gap",
            "qmp_minus_clr_z_gap_lcpm",
            "qmp_minus_clr_z_gap_metacardis",
        ),
    ]:
        pearson = stats.pearsonr(shared[lcpm_column], shared[meta_column])
        spearman = stats.spearmanr(shared[lcpm_column], shared[meta_column])
        summary_rows.append(
            {
                "comparison": name,
                "shared_exact_species": len(shared),
                "pearson_r": float(pearson.statistic),
                "spearman_rho": float(spearman.statistic),
                "interpretation": (
                    "descriptive representation sensitivity; taxa are correlated, "
                    "so no inferential P value is assigned"
                ),
            }
        )

    ordered = [
        "species",
        "feature",
        "feature_id",
        "matrix_column",
        "QMP_lcpm",
        "Row-closed_lcpm",
        "CLR_lcpm",
        "QMP_metacardis",
        "Row-closed_metacardis",
        "CLR_metacardis",
        "QMP_z_lcpm",
        "Row-closed_z_lcpm",
        "CLR_z_lcpm",
        "QMP_z_metacardis",
        "Row-closed_z_metacardis",
        "CLR_z_metacardis",
        "qmp_minus_row_closed_z_gap_lcpm",
        "qmp_minus_row_closed_z_gap_metacardis",
        "qmp_minus_clr_z_gap_lcpm",
        "qmp_minus_clr_z_gap_metacardis",
    ]
    return shared[ordered].sort_values("species"), pd.DataFrame(summary_rows)


def run(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Path]:
    association_dir = output_dir / "associations"
    prediction_dir = output_dir / "prediction"
    synthesis_dir = ensure_dir(output_dir / "synthesis")

    lcpm_all = pd.read_csv(association_dir / "lcpm_crc_vs_ctl_associations.csv")
    lcpm = lcpm_all.loc[
        lcpm_all["primary_analysis"].astype(bool)
        & lcpm_all["component"].isin(
            ["qmp", "row_closed", "clr_minimum_positive"]
        )
    ].copy()
    meta = pd.read_csv(
        association_dir / "metacardis_ihd_vs_mmc_hurdle_qmp_row_closed_clr.csv"
    )
    repeats = pd.read_csv(prediction_dir / "cv_repeat_metrics.csv")
    cv_summary = pd.read_csv(prediction_dir / "cv_summary.csv")

    microbiome_repeats = repeats.loc[
        repeats["model"].isin(PRIMARY_REPRESENTATIONS)
    ]
    microbiome_cv = cv_summary.loc[
        cv_summary["model"].isin(PRIMARY_REPRESENTATIONS)
    ].copy()
    cv_difference = pd.concat(
        [
            paired_repeat_descriptions(
                microbiome_repeats,
                cohort,
                PRIMARY_REPRESENTATIONS,
            )
            for cohort in ["LCPM", "MetaCardis"]
        ],
        ignore_index=True,
    )
    clr_cv = cv_summary.loc[
        cv_summary["model"].isin(["CLR", "CLR (multiplicative)"])
    ].copy()
    clr_difference = pd.concat(
        [
            paired_repeat_descriptions(
                repeats,
                cohort,
                ["CLR", "CLR (multiplicative)"],
                metrics=["roc_auc"],
            )
            for cohort in ["LCPM", "MetaCardis"]
        ],
        ignore_index=True,
    )

    for key, expected in EXPECTED_AUCS.items():
        observed = float(
            microbiome_cv.loc[
                (microbiome_cv["cohort"] == key[0])
                & (microbiome_cv["model"] == key[1]),
                "roc_auc_mean",
            ].iloc[0]
        )
        require(
            abs(observed - expected) < 5e-4,
            f"Unexpected AUC for {key}: expected about {expected}, observed {observed}",
        )

    da_rows = []
    for (specification, component), group in lcpm_all.groupby(
        ["filter_specification", "component"], sort=False
    ):
        da_rows.append(
            {
                "cohort": "LCPM",
                "contrast": "CRC vs CTL",
                "adjustment": "unadjusted",
                "filter_specification": specification,
                "component": component,
                "representation": group["representation"].iloc[0],
                "clr_zero_replacement": group["clr_zero_replacement"].iloc[0],
                "tested_features": len(group),
                "significant_q_lt_0_05": int(group["significant_q_lt_0_05"].sum()),
                "effect_definition": "rank-biserial CRC minus CTL",
            }
        )
    for (adjustment, component), group in meta.groupby(["adjustment", "component"]):
        da_rows.append(
            {
                "cohort": "MetaCardis",
                "contrast": "IHD372 vs MMC372",
                "adjustment": adjustment,
                "filter_specification": "pooled_outcome_blind",
                "component": component,
                "representation": _meta_representation(component),
                "clr_zero_replacement": (
                    "minimum_positive"
                    if component == "clr_minimum_positive"
                    else "multiplicative_delta_1_over_p_squared"
                    if component == "clr_multiplicative"
                    else "not_applicable"
                ),
                "tested_features": len(group),
                "significant_q_lt_0_05": int(group["significant_q_lt_0_05"].sum()),
                "effect_definition": group["effect_definition"].iloc[0],
            }
        )
    da_counts = pd.DataFrame(da_rows)

    lcpm_filter_rows = []
    for component in [
        "qmp",
        "row_closed",
        "clr_minimum_positive",
        "clr_multiplicative",
    ]:
        primary_group = lcpm_all.loc[
            (lcpm_all["filter_specification"] == "pooled_outcome_blind")
            & (lcpm_all["component"] == component)
        ]
        source_group = lcpm_all.loc[
            (lcpm_all["filter_specification"] == "source_aligned_group_union")
            & (lcpm_all["component"] == component)
        ]
        primary_calls = set(
            primary_group.loc[primary_group["significant_q_lt_0_05"], "feature"]
        )
        source_calls = set(
            source_group.loc[source_group["significant_q_lt_0_05"], "feature"]
        )
        lcpm_filter_rows.append(
            {
                "component": component,
                "pooled_outcome_blind_tested": len(primary_group),
                "source_aligned_tested": len(source_group),
                "tested_feature_overlap": len(
                    set(primary_group["feature"]) & set(source_group["feature"])
                ),
                "pooled_outcome_blind_calls": len(primary_calls),
                "source_aligned_calls": len(source_calls),
                "significant_call_overlap": len(primary_calls & source_calls),
            }
        )
    lcpm_filter_sensitivity = pd.DataFrame(lcpm_filter_rows)

    stability_rows = _stability_rows(
        "LCPM",
        "CRC vs CTL",
        "unadjusted; pooled outcome-blind filter",
        lcpm,
        "feature",
        "representation",
    )
    for adjustment in ["core", "core_plus_medications"]:
        current = meta.loc[
            (meta["adjustment"] == adjustment)
            & meta["component"].isin(
                ["qmp_nonzero", "row_closed_nonzero", "clr_minimum_positive"]
            )
        ].copy()
        current["representation"] = current["component"].map(_meta_representation)
        stability_rows.extend(
            _stability_rows(
                "MetaCardis",
                "IHD372 vs MMC372",
                adjustment,
                current,
                "matrix_column",
                "representation",
            )
        )
    stability = pd.DataFrame(stability_rows)
    shared, shared_summary = _shared_species(lcpm, meta)

    paths = {
        "cv": synthesis_dir / "microbiome_cv_summary.csv",
        "cv_differences": synthesis_dir
        / "microbiome_cv_paired_descriptive_differences.csv",
        "clr_cv": synthesis_dir / "clr_zero_replacement_cv_summary.csv",
        "clr_cv_differences": synthesis_dir
        / "clr_zero_replacement_paired_descriptive_differences.csv",
        "da_counts": synthesis_dir / "differential_association_counts.csv",
        "lcpm_filter_sensitivity": synthesis_dir
        / "lcpm_outcome_blind_filter_sensitivity.csv",
        "effect_stability": synthesis_dir / "effect_stability.csv",
        "shared_species": synthesis_dir / "shared_exact_species.csv",
        "shared_summary": synthesis_dir / "shared_species_gap_summary.csv",
        "headline": synthesis_dir / "headline_results.json",
    }
    microbiome_cv.to_csv(paths["cv"], index=False)
    cv_difference.to_csv(paths["cv_differences"], index=False)
    clr_cv.to_csv(paths["clr_cv"], index=False)
    clr_difference.to_csv(paths["clr_cv_differences"], index=False)
    da_counts.to_csv(paths["da_counts"], index=False)
    lcpm_filter_sensitivity.to_csv(paths["lcpm_filter_sensitivity"], index=False)
    stability.to_csv(paths["effect_stability"], index=False)
    shared.to_csv(paths["shared_species"], index=False)
    shared_summary.to_csv(paths["shared_summary"], index=False)

    aucs = {
        f"{cohort}__{model}": float(value)
        for cohort, model, value in microbiome_cv[
            ["cohort", "model", "roc_auc_mean"]
        ].itertuples(index=False, name=None)
    }
    qmp_row_closed_stability = stability.loc[
        (stability["representation_a"] == "QMP")
        & (stability["representation_b"] == "Row-closed")
        & stability["adjustment"].isin(
            ["unadjusted; pooled outcome-blind filter", "core"]
        )
    ]
    write_json(
        {
            "microbiome_auc_means": aucs,
            "clr_multiplicative_auc_means": {
                row.cohort: float(row.roc_auc_mean)
                for row in clr_cv.loc[
                    clr_cv["model"] == "CLR (multiplicative)"
                ].itertuples()
            },
            "shared_exact_species": len(shared),
            "qmp_row_closed_effect_correlations": {
                row.cohort: float(row.effect_pearson_r)
                for row in qmp_row_closed_stability.itertuples()
            },
            "interpretation_boundary": (
                "Do not pool diseases, raw QMP scales, or disease effects. "
                "Shared-species results compare standardized within-cohort "
                "representation sensitivity only."
            ),
        },
        paths["headline"],
    )
    print(f"synthesis complete: {len(shared)} one-to-one shared species")
    return paths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    run(args.output_dir)


if __name__ == "__main__":
    main()
