"""Part 06 — leakage-controlled repeated-CV prediction models.

The scientific comparison is deliberately paired: every representation uses
the same participant split and the same feature set within a fold.  Taxon
prevalence filtering and CLR zero replacement are learned from the training
fold only.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.utils.class_weight import compute_sample_weight

from .common import (
    classification_metrics,
    close_rows,
    clr_transform,
    ensure_dir,
    multiplicative_clr,
    paired_repeat_descriptions,
    prevalence_filter,
    read_matrix,
    require,
    summarize_repeat_metrics,
)
from .config import (
    DEFAULT_DATA_DIR,
    DEFAULT_OUTPUT_DIR,
    HGB_PARAMETERS,
    LCPM_DETECTION_LIMIT,
    N_REPEATS,
    N_SPLITS,
    PRIMARY_PREVALENCE,
    SEED,
)


MICROBIOME_MODELS = ["QMP", "Row-closed", "CLR"]
CLR_SENSITIVITY_MODEL = "CLR (multiplicative)"
CLINICAL_COLUMNS = [
    "age",
    "bmi",
    "bmi_missing",
    "male",
    "france",
    "diabetes",
    "antidiabetic_drug",
    "antihypertensive_drug",
    "lipid_lowering_drug",
    "proton_pump_inhibitor",
]


def _hgb_probability(
    train_x: pd.DataFrame,
    train_y: pd.Series,
    test_x: pd.DataFrame,
) -> np.ndarray:
    model = HistGradientBoostingClassifier(**HGB_PARAMETERS)
    weights = compute_sample_weight(class_weight="balanced", y=train_y)
    model.fit(train_x, train_y, sample_weight=weights)
    return model.predict_proba(test_x)[:, 1]


def _clinical_design(metadata: pd.DataFrame, ids: pd.Index) -> pd.DataFrame:
    """Create the prespecified MetaCardis clinical design.

    BMI is left missing here so its imputation value can be learned inside each
    training fold.  A missingness indicator remains in the model, and
    nationality is binary in this two-country primary cohort (France versus
    Denmark).
    """
    source = metadata.loc[ids]
    design = pd.DataFrame(index=ids)
    design["age"] = pd.to_numeric(source["age"], errors="raise")
    bmi = pd.to_numeric(source["bmi"], errors="coerce")
    design["bmi_missing"] = bmi.isna().astype(float)
    design["bmi"] = bmi
    design["male"] = pd.to_numeric(source["male"], errors="raise")
    design["france"] = source["nationality"].eq("France").astype(float)
    design["diabetes"] = pd.to_numeric(source["diabetes"], errors="raise")
    for column in [
        "antidiabetic_drug",
        "antihypertensive_drug",
        "lipid_lowering_drug",
        "proton_pump_inhibitor",
    ]:
        design[column] = pd.to_numeric(source[column], errors="raise")
    return design[CLINICAL_COLUMNS]


def _clinical_probability(
    train_x: pd.DataFrame,
    train_y: pd.Series,
    test_x: pd.DataFrame,
) -> np.ndarray:
    model = LogisticRegression(
        class_weight="balanced",
        max_iter=5000,
        random_state=SEED,
    )
    model.fit(train_x, train_y)
    return model.predict_proba(test_x)[:, 1]


def _fold_metric_row(
    cohort: str,
    repeat: int,
    fold: int,
    model: str,
    y_true: pd.Series,
    probability: np.ndarray,
) -> dict[str, object]:
    prediction = probability >= 0.5
    return {
        "cohort": cohort,
        "repeat": repeat,
        "fold": fold,
        "model": model,
        "n_test": len(y_true),
        "roc_auc": float(roc_auc_score(y_true, probability)),
        "average_precision": float(average_precision_score(y_true, probability)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, prediction)),
        "brier_score": float(brier_score_loss(y_true, probability)),
    }


def _run_cohort(
    cohort: str,
    metadata: pd.DataFrame,
    qmp: pd.DataFrame,
    candidate_features: list[str],
    ids: pd.Index,
    labels: pd.Series,
    detection_limit: float,
    include_clinical: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Preserve the published matrix order.  The random splitter acts on row
    # positions, so this is part of exact reproducibility.
    ids = qmp.index[qmp.index.isin(ids)]
    qmp = qmp.loc[ids]
    labels = labels.loc[ids].astype(int)
    row_closed = close_rows(qmp)
    filter_reference = row_closed if cohort == "LCPM" else qmp
    clinical = _clinical_design(metadata, ids) if include_clinical else None

    model_names = [*MICROBIOME_MODELS, CLR_SENSITIVITY_MODEL]
    if include_clinical:
        model_names.extend(
            [
                "Clinical (logistic)",
                "Clinical (HGB)",
                "QMP + clinical",
                "Row-closed + clinical",
                "CLR + clinical",
            ]
        )
    probabilities = {
        model: np.full((N_REPEATS, len(ids)), np.nan, dtype=float)
        for model in model_names
    }
    fold_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    selected = Counter()

    splitter = RepeatedStratifiedKFold(
        n_splits=N_SPLITS,
        n_repeats=N_REPEATS,
        random_state=SEED,
    )
    for iteration, (train_position, test_position) in enumerate(
        splitter.split(qmp, labels)
    ):
        repeat = iteration // N_SPLITS + 1
        fold = iteration % N_SPLITS + 1
        train_ids = ids[train_position]
        test_ids = ids[test_position]
        train_y = labels.iloc[train_position]
        test_y = labels.iloc[test_position]
        features = prevalence_filter(
            filter_reference.loc[train_ids, candidate_features],
            train_y,
            PRIMARY_PREVALENCE,
            detection_limit,
        )
        require(bool(features), f"{cohort} repeat {repeat} fold {fold}: no features")
        selected.update(features)

        clr_train, clr_test, replacements = clr_transform(
            row_closed.loc[train_ids, features],
            row_closed.loc[test_ids, features],
        )
        multiplicative_clr_train = multiplicative_clr(
            row_closed.loc[train_ids, features]
        )
        multiplicative_clr_test = multiplicative_clr(
            row_closed.loc[test_ids, features]
        )
        microbiome = {
            "QMP": (qmp.loc[train_ids, features], qmp.loc[test_ids, features]),
            "Row-closed": (
                row_closed.loc[train_ids, features],
                row_closed.loc[test_ids, features],
            ),
            "CLR": (clr_train, clr_test),
            CLR_SENSITIVITY_MODEL: (
                multiplicative_clr_train,
                multiplicative_clr_test,
            ),
        }

        fold_predictions: dict[str, np.ndarray] = {}
        for model_name, (train_x, test_x) in microbiome.items():
            fold_predictions[model_name] = _hgb_probability(train_x, train_y, test_x)

        clinical_bmi_median = np.nan
        if clinical is not None:
            clinical_train = clinical.loc[train_ids].copy()
            clinical_test = clinical.loc[test_ids].copy()
            clinical_bmi_median = float(clinical_train["bmi"].median())
            clinical_train["bmi"] = clinical_train["bmi"].fillna(
                clinical_bmi_median
            )
            clinical_test["bmi"] = clinical_test["bmi"].fillna(
                clinical_bmi_median
            )
            fold_predictions["Clinical (logistic)"] = _clinical_probability(
                clinical_train, train_y, clinical_test
            )
            fold_predictions["Clinical (HGB)"] = _hgb_probability(
                clinical_train, train_y, clinical_test
            )
            for model_name in MICROBIOME_MODELS:
                train_x, test_x = microbiome[model_name]
                combined_train = train_x.join(
                    clinical_train.add_prefix("clinical__")
                )
                combined_test = test_x.join(clinical_test.add_prefix("clinical__"))
                fold_predictions[f"{model_name} + clinical"] = _hgb_probability(
                    combined_train, train_y, combined_test
                )

        for model_name, probability in fold_predictions.items():
            probabilities[model_name][repeat - 1, test_position] = probability
            fold_rows.append(
                _fold_metric_row(
                    cohort, repeat, fold, model_name, test_y, probability
                )
            )

        audit_rows.append(
            {
                "cohort": cohort,
                "repeat": repeat,
                "fold": fold,
                "training_samples": len(train_ids),
                "test_samples": len(test_ids),
                "training_positive": int(train_y.sum()),
                "training_negative": int((1 - train_y).sum()),
                "selected_features": len(features),
                "minimum_clr_replacement": float(replacements.min()),
                "maximum_clr_replacement": float(replacements.max()),
                "multiplicative_clr_delta": float(1.0 / (len(features) ** 2)),
                "clinical_training_bmi_median": clinical_bmi_median,
            }
        )
        if (iteration + 1) % 10 == 0:
            print(f"{cohort}: completed {iteration + 1}/{N_SPLITS * N_REPEATS} folds", flush=True)

    repeat_rows: list[dict[str, object]] = []
    for model_name, matrix in probabilities.items():
        require(np.isfinite(matrix).all(), f"Missing OOF predictions for {cohort} {model_name}")
        for repeat_position in range(N_REPEATS):
            row: dict[str, object] = {
                "cohort": cohort,
                "model": model_name,
                "repeat": repeat_position + 1,
                "folds": N_SPLITS,
                "n": len(labels),
                "n_negative": int((1 - labels).sum()),
                "n_positive": int(labels.sum()),
            }
            row.update(classification_metrics(labels.to_numpy(), matrix[repeat_position]))
            repeat_rows.append(row)

    selection_rows = [
        {
            "cohort": cohort,
            "feature": feature,
            "selected_folds": selected.get(feature, 0),
            "total_folds": N_SPLITS * N_REPEATS,
            "selection_frequency": selected.get(feature, 0) / (N_SPLITS * N_REPEATS),
        }
        for feature in candidate_features
    ]
    return (
        pd.DataFrame(repeat_rows),
        pd.DataFrame(fold_rows),
        pd.DataFrame(audit_rows),
        pd.DataFrame(selection_rows),
    )


def run(
    data_dir: Path = DEFAULT_DATA_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    include_metacardis_clinical: bool = True,
) -> dict[str, Path]:
    result_dir = ensure_dir(output_dir / "prediction")

    lcpm_dir = data_dir / "processed" / "lcpm"
    lcpm_metadata = pd.read_csv(lcpm_dir / "metadata.csv.gz").set_index(
        "participant_id"
    )
    lcpm_taxonomy = pd.read_csv(lcpm_dir / "taxonomy.csv.gz")
    lcpm_qmp = read_matrix(lcpm_dir / "qmp_cells_per_g.csv")
    lcpm_ids = lcpm_metadata.index[
        lcpm_metadata["diagnosis"].isin(["CTL", "CRC"])
    ]
    lcpm_labels = lcpm_metadata.loc[lcpm_ids, "diagnosis"].eq("CRC").astype(int)
    lcpm_candidates = lcpm_taxonomy.loc[
        lcpm_taxonomy["clean_candidate"], "feature_id"
    ].tolist()
    lcpm_results = _run_cohort(
        "LCPM",
        lcpm_metadata,
        lcpm_qmp,
        lcpm_candidates,
        lcpm_ids,
        lcpm_labels,
        LCPM_DETECTION_LIMIT,
        include_clinical=False,
    )

    meta_dir = data_dir / "processed" / "metacardis"
    meta_metadata = pd.read_csv(meta_dir / "metadata.csv.gz").set_index(
        "participant_id"
    )
    meta_taxonomy = pd.read_csv(meta_dir / "taxonomy.csv.gz")
    meta_qmp = read_matrix(meta_dir / "qmp_index.csv")
    meta_mask = meta_metadata["quantitative_profile_available"].astype(bool) & (
        meta_metadata["ihd_member"].astype(bool)
        | meta_metadata["mmc_member"].astype(bool)
    )
    meta_ids = meta_metadata.index[meta_mask]
    meta_labels = meta_metadata.loc[meta_ids, "ihd_member"].astype(int)
    meta_candidates = meta_taxonomy.loc[
        meta_taxonomy["clean_candidate"], "matrix_column"
    ].tolist()
    meta_results = _run_cohort(
        "MetaCardis",
        meta_metadata,
        meta_qmp,
        meta_candidates,
        meta_ids,
        meta_labels,
        np.nextafter(0.0, 1.0),
        include_clinical=include_metacardis_clinical,
    )

    repeats = pd.concat([lcpm_results[0], meta_results[0]], ignore_index=True)
    fold_metrics = pd.concat([lcpm_results[1], meta_results[1]], ignore_index=True)
    fold_audit = pd.concat([lcpm_results[2], meta_results[2]], ignore_index=True)
    selection = pd.concat([lcpm_results[3], meta_results[3]], ignore_index=True)
    summary = summarize_repeat_metrics(repeats, ["cohort", "model"])

    comparison_frames = []
    for cohort, cohort_models in repeats.groupby("cohort", sort=False)["model"]:
        comparison_frames.append(
            paired_repeat_descriptions(
                repeats,
                cohort,
                cohort_models.drop_duplicates().tolist(),
            )
        )
    comparisons = pd.concat(comparison_frames, ignore_index=True)

    paths = {
        "repeats": result_dir / "cv_repeat_metrics.csv",
        "summary": result_dir / "cv_summary.csv",
        "comparisons": result_dir / "cv_paired_descriptive_differences.csv",
        "fold_metrics": result_dir / "cv_fold_metrics.csv",
        "fold_audit": result_dir / "cv_fold_audit.csv",
        "selection": result_dir / "cv_feature_selection_frequency.csv",
    }
    repeats.to_csv(paths["repeats"], index=False)
    summary.to_csv(paths["summary"], index=False)
    comparisons.to_csv(paths["comparisons"], index=False)
    fold_metrics.to_csv(paths["fold_metrics"], index=False)
    fold_audit.to_csv(paths["fold_audit"], index=False)
    selection.to_csv(paths["selection"], index=False)

    headline = summary.loc[
        summary["model"].isin(MICROBIOME_MODELS),
        ["cohort", "model", "roc_auc_mean"],
    ]
    print("Primary repeated-CV AUCs:")
    print(headline.to_string(index=False))
    return paths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--microbiome-only",
        action="store_true",
        help="Skip the optional MetaCardis clinical and combined models.",
    )
    args = parser.parse_args()
    run(
        args.data_dir,
        args.output_dir,
        include_metacardis_clinical=not args.microbiome_only,
    )


if __name__ == "__main__":
    main()
