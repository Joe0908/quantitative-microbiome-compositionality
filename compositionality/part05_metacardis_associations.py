"""Part 05 — MetaCardis prevalence, abundance, and CLR association models."""

from __future__ import annotations

import argparse
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from .common import (
    bh_fdr,
    close_rows,
    ensure_dir,
    fixed_clr,
    multiplicative_clr,
    pooled_prevalence_filter,
    read_matrix,
    require,
    write_json,
)
from .config import DEFAULT_DATA_DIR, DEFAULT_OUTPUT_DIR, PRIMARY_PREVALENCE


def _base_covariates(metadata: pd.DataFrame, disease: pd.Series) -> pd.DataFrame:
    result = pd.DataFrame(index=metadata.index)
    result["disease"] = disease.astype(float)
    result["age"] = pd.to_numeric(metadata["age"], errors="coerce")
    bmi = pd.to_numeric(metadata["bmi"], errors="coerce")
    result["bmi_missing"] = bmi.isna().astype(float)
    result["bmi"] = bmi.fillna(bmi.median())
    result["male"] = pd.to_numeric(metadata["male"], errors="coerce")
    result["diabetes"] = pd.to_numeric(metadata["diabetes"], errors="coerce")
    nationality = pd.get_dummies(
        metadata["nationality"].astype("category"),
        prefix="nationality",
        drop_first=True,
        dtype=float,
    )
    result = result.join(nationality)
    for column in [
        "antidiabetic_drug",
        "antihypertensive_drug",
        "lipid_lowering_drug",
        "proton_pump_inhibitor",
    ]:
        result[column] = pd.to_numeric(metadata[column], errors="coerce")
    return result


def _design(covariates: pd.DataFrame, index: pd.Index, nuisance: list[str]) -> pd.DataFrame | None:
    subset = covariates.loc[index]
    if subset["disease"].nunique() < 2:
        return None
    design = pd.DataFrame({"intercept": 1.0, "disease": subset["disease"]}, index=index)
    current_rank = np.linalg.matrix_rank(design.to_numpy(dtype=float))
    for column in nuisance:
        candidate = subset[column].astype(float)
        if candidate.isna().any() or candidate.nunique() <= 1:
            continue
        proposed = design.assign(**{column: candidate})
        proposed_rank = np.linalg.matrix_rank(proposed.to_numpy(dtype=float))
        if proposed_rank > current_rank:
            design = proposed
            current_rank = proposed_rank
    return design


def _fit_logistic(response: pd.Series, design: pd.DataFrame | None) -> tuple[float, float, float, str]:
    if design is None or response.nunique() < 2:
        return math.nan, math.nan, math.nan, "not_estimable"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = sm.GLM(response.astype(float), design, family=sm.families.Binomial()).fit(
                maxiter=200, disp=0
            )
        values = (fit.params["disease"], fit.bse["disease"], fit.pvalues["disease"])
        if not np.isfinite(values).all():
            return math.nan, math.nan, math.nan, "non_finite"
        return float(values[0]), float(values[1]), float(values[2]), "ok"
    except Exception as error:  # explicit status is retained in the result table
        return math.nan, math.nan, math.nan, f"failed:{type(error).__name__}"


def _fit_ols(response: pd.Series, design: pd.DataFrame | None) -> tuple[float, float, float, str]:
    if design is None:
        return math.nan, math.nan, math.nan, "not_estimable"
    try:
        # ``use_t=True`` reports finite-sample t tests with the residual
        # degrees of freedom.  This is the specification used in the
        # completed analysis (the HC3 coefficients and standard errors are
        # unchanged, but the p values are slightly less anti-conservative
        # than the asymptotic normal approximation).
        columns = list(design.columns)
        nuisance = [column for column in columns if column not in {"intercept", "disease"}]
        # HC3 contains (1 - leverage)^-2.  A single participant carrying a
        # sparse missingness indicator can have leverage numerically equal to
        # one, and the last bit depends on otherwise irrelevant column order.
        # Retry algebraically identical orders so a 0/0 rounding artefact does
        # not erase a valid coefficient and test.
        orders = [
            columns,
            ["intercept", "disease", *reversed(nuisance)],
            ["intercept", *nuisance, "disease"],
            ["disease", "intercept", *nuisance],
            [*nuisance, "intercept", "disease"],
            [*nuisance, "disease", "intercept"],
            list(reversed(columns)),
        ]
        for offset in range(len(nuisance)):
            orders.append(
                [
                    "intercept",
                    "disease",
                    *nuisance[offset:],
                    *nuisance[:offset],
                ]
            )
        for order in orders:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                fit = sm.OLS(response.astype(float), design[order]).fit(
                    cov_type="HC3", use_t=True
                )
            values = (
                fit.params["disease"],
                fit.bse["disease"],
                fit.pvalues["disease"],
            )
            if np.isfinite(values).all():
                return float(values[0]), float(values[1]), float(values[2]), "ok"
        return math.nan, math.nan, math.nan, "non_finite"
    except Exception as error:
        return math.nan, math.nan, math.nan, f"failed:{type(error).__name__}"


def run(
    data_dir: Path = DEFAULT_DATA_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> dict[str, Path]:
    processed = data_dir / "processed" / "metacardis"
    results_dir = ensure_dir(output_dir / "associations")
    metadata = pd.read_csv(processed / "metadata.csv.gz").set_index("participant_id")
    taxonomy = pd.read_csv(processed / "taxonomy.csv.gz")
    qmp = read_matrix(processed / "qmp_index.csv")
    rmp = close_rows(qmp)

    primary_mask = metadata["quantitative_profile_available"].astype(bool) & (
        metadata["ihd_member"].astype(bool) | metadata["mmc_member"].astype(bool)
    )
    primary_metadata = metadata.loc[primary_mask].copy()
    disease = primary_metadata["ihd_member"].astype(int)
    require(len(primary_metadata) == 672, "Expected 672 MetaCardis primary profiles")
    require(int(disease.sum()) == 303, "Expected 303 IHD profiles")

    candidate_columns = taxonomy.loc[taxonomy["clean_candidate"], "matrix_column"].tolist()
    qmp = qmp.loc[primary_metadata.index, candidate_columns]
    rmp = rmp.loc[primary_metadata.index, candidate_columns]
    # Pooled prevalence is independent of IHD/MMC labels and therefore does
    # not select hypotheses using the outcome later tested by the models.
    features = pooled_prevalence_filter(
        qmp,
        PRIMARY_PREVALENCE,
        np.nextafter(0.0, 1.0),
    )
    require(
        len(features) == 404,
        f"Expected 404 pooled outcome-blind MetaCardis features, observed {len(features)}",
    )
    qmp = qmp[features]
    rmp = rmp[features]
    clr_minimum = fixed_clr(rmp)
    clr_multiplicative = multiplicative_clr(rmp)

    feature_info = taxonomy.set_index("matrix_column").loc[features]
    covariates = _base_covariates(primary_metadata, disease)
    nationality_columns = [column for column in covariates if column.startswith("nationality_")]
    # Keep the prespecified nuisance order.  Besides making the model record
    # readable, the order avoids machine-precision HC3 leverage artefacts for a
    # handful of sparse positive-only fits.
    core = [*nationality_columns, "age", "bmi", "bmi_missing", "male", "diabetes"]
    medications = [
        "antidiabetic_drug",
        "antihypertensive_drug",
        "lipid_lowering_drug",
        "proton_pump_inhibitor",
    ]
    adjustments = {
        "unadjusted": [],
        "core": core,
        "core_plus_medications": [*core, *medications],
    }

    rows: list[dict[str, object]] = []
    for adjustment, nuisance in adjustments.items():
        full_design = _design(covariates, primary_metadata.index, nuisance)
        for position, feature in enumerate(features, start=1):
            presence = qmp[feature] > 0
            positive_index = presence.index[presence]
            positive_design = _design(covariates, positive_index, nuisance)
            component_inputs = [
                ("prevalence", presence.astype(float), full_design, _fit_logistic, "log odds IHD vs MMC"),
                ("qmp_nonzero", np.log(qmp.loc[positive_index, feature]), positive_design, _fit_ols, "log non-zero QMP IHD minus MMC"),
                ("row_closed_nonzero", np.log(rmp.loc[positive_index, feature]), positive_design, _fit_ols, "log non-zero QMP-derived row-closed abundance IHD minus MMC"),
                ("clr_minimum_positive", clr_minimum[feature], full_design, _fit_ols, "CLR IHD minus MMC; minimum-positive zero replacement"),
                ("clr_multiplicative", clr_multiplicative[feature], full_design, _fit_ols, "CLR IHD minus MMC; multiplicative zero replacement"),
            ]
            for component, response, design, fit_function, definition in component_inputs:
                effect, standard_error, p_value, status = fit_function(response, design)
                rows.append(
                    {
                        "cohort": "MetaCardis",
                        "contrast": "IHD372 vs MMC372",
                        "adjustment": adjustment,
                        "component": component,
                        "matrix_column": feature,
                        "feature_id": feature_info.loc[feature, "Feature ID"],
                        "species": feature_info.loc[feature, "species"],
                        "effect": effect,
                        "standard_error": standard_error,
                        "p_value": p_value,
                        "effect_definition": definition,
                        "n_total": len(primary_metadata),
                        "n_analyzed": len(response),
                        "ihd_present": int(presence.loc[disease.index[disease == 1]].sum()),
                        "mmc_present": int(presence.loc[disease.index[disease == 0]].sum()),
                        "status": status,
                    }
                )
            if position % 50 == 0:
                print(f"{adjustment}: fitted {position}/{len(features)} taxa", flush=True)

    results = pd.DataFrame(rows)
    results["q_value_bh"] = np.nan
    for _, index in results.groupby(["adjustment", "component"]).groups.items():
        results.loc[index, "q_value_bh"] = bh_fdr(results.loc[index, "p_value"])
    results["significant_q_lt_0_05"] = results["q_value_bh"] < 0.05

    result_path = results_dir / "metacardis_ihd_vs_mmc_hurdle_qmp_row_closed_clr.csv"
    summary_path = results_dir / "metacardis_association_summary.json"
    results.to_csv(result_path, index=False)
    counts = (
        results.groupby(["adjustment", "component"])["significant_q_lt_0_05"]
        .sum()
        .astype(int)
        .unstack(fill_value=0)
    )
    core_results = results.loc[results["adjustment"] == "core"]
    qmp_effect = core_results.loc[core_results["component"] == "qmp_nonzero"].set_index("matrix_column")["effect"]
    row_closed_effect = core_results.loc[
        core_results["component"] == "row_closed_nonzero"
    ].set_index("matrix_column")["effect"]
    effect_correlation = float(qmp_effect.corr(row_closed_effect))
    direction_agreement = float(
        (np.sign(qmp_effect) == np.sign(row_closed_effect)).mean()
    )
    write_json(
        {
            "contrast": "IHD372 vs MMC372",
            "samples": 672,
            "ihd": 303,
            "mmc": 369,
            "filter": "pooled outcome-blind prevalence >=5%",
            "tested_features": 404,
            "bmi_median_imputation": float(pd.to_numeric(primary_metadata["bmi"], errors="coerce").median()),
            "significant_calls": counts.to_dict(orient="index"),
            "core_qmp_row_closed_effect_correlation": effect_correlation,
            "core_qmp_row_closed_direction_agreement": direction_agreement,
            "medication_adjustment_interpretation": (
                "sensitivity analysis under an alternative covariate set; "
                "not a causal estimate of medication effects"
            ),
        },
        summary_path,
    )
    print("MetaCardis association calls:")
    print(counts.to_string())
    return {"results": result_path, "summary": summary_path}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    run(args.data_dir, args.output_dir)


if __name__ == "__main__":
    main()
