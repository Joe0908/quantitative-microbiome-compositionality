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
    paired_repeat_comparisons,
    require,
    write_json,
    zscore,
)
from .config import DEFAULT_OUTPUT_DIR


EXPECTED_AUCS = {
    ("Galazzo/LCPM", "QMP"): 0.6589828749351323,
    ("Galazzo/LCPM", "RMP"): 0.6528697457187339,
    ("Galazzo/LCPM", "CLR"): 0.6416917488323820,
    ("MetaCardis", "QMP"): 0.6393642616294150,
    ("MetaCardis", "RMP"): 0.6472027690573936,
    ("MetaCardis", "CLR"): 0.6428577817131307,
}


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
        for representation in ["QMP", "RMP", "CLR"]
    }
    rows = []
    for left, right in combinations(["QMP", "RMP", "CLR"], 2):
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
        "rmp_nonzero": "RMP",
        "clr": "CLR",
    }.get(value)


def _shared_species(
    lcpm: pd.DataFrame, meta: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    lcpm_effect = lcpm.pivot(
        index="feature", columns="representation", values="effect"
    )[["QMP", "RMP", "CLR"]]
    for representation in ["QMP", "RMP", "CLR"]:
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
        & meta["component"].isin(["qmp_nonzero", "rmp_nonzero", "clr"])
    ].copy()
    meta_core["representation"] = meta_core["component"].map(_meta_representation)
    meta_effect = meta_core.pivot(
        index=["matrix_column", "feature_id", "species"],
        columns="representation",
        values="effect",
    ).reset_index()
    for representation in ["QMP", "RMP", "CLR"]:
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

    require(len(lcpm_one) == 109, f"Expected 109 one-to-one LCPM species; got {len(lcpm_one)}")
    require(len(meta_strict) == 145, f"Expected 145 strict MetaCardis rows; got {len(meta_strict)}")
    require(len(duplicated_labels) == 14, f"Expected 14 duplicated MetaCardis species; got {len(duplicated_labels)}")
    require(int(meta_counts.loc[duplicated_labels].sum()) == 32, "Expected 32 duplicated MetaCardis MGS rows")
    require(len(meta_one) == 113, f"Expected 113 one-to-one MetaCardis species; got {len(meta_one)}")

    shared = lcpm_one.merge(
        meta_one,
        on="species",
        how="inner",
        suffixes=("_lcpm", "_metacardis"),
        validate="one_to_one",
    )
    require(len(shared) == 51, f"Expected 51 shared exact species; got {len(shared)}")
    for cohort_suffix in ["lcpm", "metacardis"]:
        shared[f"qmp_minus_rmp_z_gap_{cohort_suffix}"] = (
            shared[f"QMP_z_{cohort_suffix}"] - shared[f"RMP_z_{cohort_suffix}"]
        )
        shared[f"qmp_minus_clr_z_gap_{cohort_suffix}"] = (
            shared[f"QMP_z_{cohort_suffix}"] - shared[f"CLR_z_{cohort_suffix}"]
        )

    summary_rows = []
    for name, lcpm_column, meta_column in [
        (
            "QMP minus RMP standardized gap",
            "qmp_minus_rmp_z_gap_lcpm",
            "qmp_minus_rmp_z_gap_metacardis",
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
                "pearson_p": float(pearson.pvalue),
                "spearman_rho": float(spearman.statistic),
                "spearman_p": float(spearman.pvalue),
                "interpretation": "representation sensitivity, not disease-effect replication",
            }
        )

    ordered = [
        "species",
        "feature",
        "feature_id",
        "matrix_column",
        "QMP_lcpm",
        "RMP_lcpm",
        "CLR_lcpm",
        "QMP_metacardis",
        "RMP_metacardis",
        "CLR_metacardis",
        "QMP_z_lcpm",
        "RMP_z_lcpm",
        "CLR_z_lcpm",
        "QMP_z_metacardis",
        "RMP_z_metacardis",
        "CLR_z_metacardis",
        "qmp_minus_rmp_z_gap_lcpm",
        "qmp_minus_rmp_z_gap_metacardis",
        "qmp_minus_clr_z_gap_lcpm",
        "qmp_minus_clr_z_gap_metacardis",
    ]
    return shared[ordered].sort_values("species"), pd.DataFrame(summary_rows)


def run(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Path]:
    association_dir = output_dir / "associations"
    prediction_dir = output_dir / "prediction"
    synthesis_dir = ensure_dir(output_dir / "synthesis")

    lcpm = pd.read_csv(association_dir / "lcpm_crc_vs_ctl_qmp_rmp_clr.csv")
    meta = pd.read_csv(
        association_dir / "metacardis_ihd_vs_mmc_hurdle_qmp_rmp_clr.csv"
    )
    repeats = pd.read_csv(prediction_dir / "cv_repeat_metrics.csv")
    cv_summary = pd.read_csv(prediction_dir / "cv_summary.csv")

    microbiome_repeats = repeats.loc[repeats["model"].isin(["QMP", "RMP", "CLR"])]
    microbiome_cv = cv_summary.loc[
        cv_summary["model"].isin(["QMP", "RMP", "CLR"])
    ].copy()
    cv_difference = pd.concat(
        [
            paired_repeat_comparisons(
                microbiome_repeats,
                cohort,
                ["QMP", "RMP", "CLR"],
            )
            for cohort in ["Galazzo/LCPM", "MetaCardis"]
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
    for representation, group in lcpm.groupby("representation"):
        da_rows.append(
            {
                "cohort": "Galazzo/LCPM",
                "contrast": "CRC vs CTL",
                "adjustment": "unadjusted pairwise",
                "component": representation,
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
                "component": component,
                "tested_features": len(group),
                "significant_q_lt_0_05": int(group["significant_q_lt_0_05"].sum()),
                "effect_definition": group["effect_definition"].iloc[0],
            }
        )
    da_counts = pd.DataFrame(da_rows)

    stability_rows = _stability_rows(
        "Galazzo/LCPM",
        "CRC vs CTL",
        "unadjusted pairwise",
        lcpm,
        "feature",
        "representation",
    )
    for adjustment in ["core", "core_plus_medications"]:
        current = meta.loc[
            (meta["adjustment"] == adjustment)
            & meta["component"].isin(["qmp_nonzero", "rmp_nonzero", "clr"])
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
        "cv_differences": synthesis_dir / "microbiome_cv_paired_differences.csv",
        "da_counts": synthesis_dir / "differential_association_counts.csv",
        "effect_stability": synthesis_dir / "effect_stability.csv",
        "shared_species": synthesis_dir / "shared_exact_species.csv",
        "shared_summary": synthesis_dir / "shared_species_gap_summary.csv",
        "headline": synthesis_dir / "headline_results.json",
    }
    microbiome_cv.to_csv(paths["cv"], index=False)
    cv_difference.to_csv(paths["cv_differences"], index=False)
    da_counts.to_csv(paths["da_counts"], index=False)
    stability.to_csv(paths["effect_stability"], index=False)
    shared.to_csv(paths["shared_species"], index=False)
    shared_summary.to_csv(paths["shared_summary"], index=False)

    aucs = {
        f"{cohort}__{model}": float(value)
        for cohort, model, value in microbiome_cv[
            ["cohort", "model", "roc_auc_mean"]
        ].itertuples(index=False, name=None)
    }
    qmp_rmp_stability = stability.loc[
        (stability["representation_a"] == "QMP")
        & (stability["representation_b"] == "RMP")
        & stability["adjustment"].isin(["unadjusted pairwise", "core"])
    ]
    write_json(
        {
            "microbiome_auc_means": aucs,
            "shared_exact_species": len(shared),
            "qmp_rmp_effect_correlations": {
                row.cohort: float(row.effect_pearson_r)
                for row in qmp_rmp_stability.itertuples()
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
