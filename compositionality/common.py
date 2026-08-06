"""Shared transformations, statistics, validation, and I/O helpers."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
)


BINOMIAL_SPACE = re.compile(r"^[A-Z][A-Za-z-]+ [a-z][a-z-]+$")
BINOMIAL_DOT = re.compile(r"^[A-Z][A-Za-z-]+\.[a-z][a-z-]+$")
CAG_SUFFIX = re.compile(r"(CAG\d+(?:_\d+)?)$")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(obj: object, path: Path) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")


def write_table(frame: pd.DataFrame, path: Path, index: bool = False) -> None:
    ensure_dir(path.parent)
    frame.to_csv(path, index=index)


def read_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    return frame.apply(pd.to_numeric, errors="raise")


def bh_fdr(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p)
    if finite.any():
        observed = p[finite]
        order = np.argsort(observed)
        ranked = observed[order]
        adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        restored = np.empty_like(adjusted)
        restored[order] = np.clip(adjusted, 0.0, 1.0)
        q[finite] = restored
    return q


def t_interval(values: Sequence[float], confidence: float = 0.95) -> tuple[float, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return (math.nan, math.nan)
    sem = stats.sem(x)
    half = stats.t.ppf((1 + confidence) / 2, x.size - 1) * sem
    return (float(x.mean() - half), float(x.mean() + half))


def close_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    row_sums = matrix.sum(axis=1)
    if (row_sums <= 0).any():
        bad = row_sums.index[row_sums <= 0].tolist()[:5]
        raise ValueError(f"Cannot row-close zero-sum profiles; examples: {bad}")
    return matrix.div(row_sums, axis=0)


def prevalence_filter(
    reference: pd.DataFrame,
    labels: pd.Series,
    threshold: float,
    detection_limit: float,
) -> list[str]:
    """Retain features detected at the requested prevalence in either class."""
    labels = labels.loc[reference.index]
    keep = np.zeros(reference.shape[1], dtype=bool)
    for group in sorted(labels.unique()):
        group_index = labels.index[labels == group]
        detected = reference.loc[group_index].to_numpy() >= detection_limit
        keep |= detected.mean(axis=0) >= threshold
    return reference.columns[keep].tolist()


def fit_clr_reference(train_rmp: pd.DataFrame) -> pd.Series:
    replacements: dict[str, float] = {}
    for column in train_rmp.columns:
        positive = train_rmp.loc[train_rmp[column] > 0, column]
        if positive.empty:
            raise ValueError(f"Feature {column!r} has no positive training values")
        replacements[column] = float(positive.min())
    return pd.Series(replacements)


def apply_clr(matrix: pd.DataFrame, replacements: pd.Series) -> pd.DataFrame:
    matrix = matrix.loc[:, replacements.index].astype(float)
    replaced = matrix.mask(matrix <= 0, replacements, axis=1)
    logged = np.log(replaced)
    return logged.sub(logged.mean(axis=1), axis=0)


def clr_transform(
    train_rmp: pd.DataFrame, test_rmp: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    replacements = fit_clr_reference(train_rmp)
    return (
        apply_clr(train_rmp, replacements),
        apply_clr(test_rmp, replacements),
        replacements,
    )


def fixed_clr(rmp: pd.DataFrame) -> pd.DataFrame:
    replacements = fit_clr_reference(rmp)
    return apply_clr(rmp, replacements)


def classification_metrics(y_true: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    prediction = (probability >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, prediction, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(y_true, probability)),
        "average_precision": float(average_precision_score(y_true, probability)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, prediction)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else math.nan,
        "specificity": float(tn / (tn + fp)) if tn + fp else math.nan,
        "brier_score": float(brier_score_loss(y_true, probability)),
    }


def summarize_repeat_metrics(repeats: pd.DataFrame, group_columns: Sequence[str]) -> pd.DataFrame:
    metric_columns = [
        "roc_auc",
        "average_precision",
        "balanced_accuracy",
        "sensitivity",
        "specificity",
        "brier_score",
    ]
    rows: list[dict[str, object]] = []
    for key, group in repeats.groupby(list(group_columns), dropna=False, sort=False):
        key = (key,) if not isinstance(key, tuple) else key
        row = dict(zip(group_columns, key))
        row["n"] = int(group["n"].iloc[0])
        row["n_negative"] = int(group["n_negative"].iloc[0])
        row["n_positive"] = int(group["n_positive"].iloc[0])
        row["cv_repeats"] = int(group["repeat"].nunique())
        row["cv_folds_per_repeat"] = int(group["folds"].iloc[0])
        for metric in metric_columns:
            values = group[metric].to_numpy(dtype=float)
            low, high = t_interval(values)
            row[f"{metric}_mean"] = float(np.mean(values))
            row[f"{metric}_sd"] = float(np.std(values, ddof=1))
            row[f"{metric}_ci95_low"] = low
            row[f"{metric}_ci95_high"] = high
        rows.append(row)
    return pd.DataFrame(rows)


def paired_repeat_comparisons(
    repeats: pd.DataFrame,
    cohort: str,
    models: Sequence[str],
    metrics: Sequence[str] = (
        "roc_auc",
        "average_precision",
        "balanced_accuracy",
        "brier_score",
    ),
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    subset = repeats.loc[repeats["cohort"] == cohort]
    for metric in metrics:
        pivot = subset.pivot(index="repeat", columns="model", values=metric)
        for i, model_a in enumerate(models):
            for model_b in models[i + 1 :]:
                difference = (pivot[model_a] - pivot[model_b]).dropna().to_numpy()
                low, high = t_interval(difference)
                if np.allclose(difference, 0):
                    pvalue = 1.0
                else:
                    pvalue = float(stats.wilcoxon(difference, alternative="two-sided").pvalue)
                rows.append(
                    {
                        "cohort": cohort,
                        "metric": metric,
                        "model_a": model_a,
                        "model_b": model_b,
                        "mean_a_minus_b": float(difference.mean()),
                        "sd_difference": float(difference.std(ddof=1)),
                        "ci95_low": low,
                        "ci95_high": high,
                        "paired_wilcoxon_p": pvalue,
                        "positive_difference_favors": (
                            model_b if metric == "brier_score" else model_a
                        ),
                    }
                )
    result = pd.DataFrame(rows)
    result["paired_wilcoxon_q_bh"] = bh_fdr(result["paired_wilcoxon_p"])
    return result


def yes_no_to_binary(series: pd.Series) -> pd.Series:
    mapping = {
        "yes": 1.0,
        "no": 0.0,
        "1": 1.0,
        "0": 0.0,
        "true": 1.0,
        "false": 0.0,
    }
    normalized = series.astype(str).str.strip().str.lower().map(mapping)
    normalized[series.isna()] = np.nan
    return normalized


def cag_key(value: object) -> str:
    match = CAG_SUFFIX.search(str(value))
    if not match:
        raise ValueError(f"No CAG suffix in {value!r}")
    return match.group(1)


def zscore(values: pd.Series) -> pd.Series:
    numeric = values.astype(float)
    return (numeric - numeric.mean()) / numeric.std(ddof=0)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def finite_or_nan(value: object) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return math.nan
    return numeric if np.isfinite(numeric) else math.nan


def ordered_unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(values))
