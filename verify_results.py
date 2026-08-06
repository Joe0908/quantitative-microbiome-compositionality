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
                cv.loc[(cv.cohort == cohort) & (cv.model == model), "roc_auc_mean"].iloc[0]
            )
            require(
                abs(observed - target) <= args.auc_tolerance,
                f"{cohort} {model} AUC: expected {target}, observed {observed}",
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

    da = pd.read_csv(args.output_dir / "synthesis" / "differential_association_counts.csv")
    lcpm_map = dict(
        da.loc[da.cohort == "Galazzo/LCPM", ["component", "significant_q_lt_0_05"]]
        .itertuples(index=False, name=None)
    )
    require(lcpm_map == expected["differential_calls_q_lt_0_05"]["Galazzo/LCPM"], "LCPM DA calls differ")

    component_names = {
        "prevalence": "prevalence",
        "qmp_nonzero": "QMP",
        "rmp_nonzero": "RMP",
        "clr": "CLR",
    }
    for adjustment, expected_key in [
        ("core", "MetaCardis_core"),
        ("core_plus_medications", "MetaCardis_core_plus_medications"),
    ]:
        current = da.loc[(da.cohort == "MetaCardis") & (da.adjustment == adjustment)].copy()
        current["label"] = current.component.map(component_names)
        observed = dict(current[["label", "significant_q_lt_0_05"]].itertuples(index=False, name=None))
        require(observed == expected["differential_calls_q_lt_0_05"][expected_key], f"{expected_key} DA calls differ")

    shared = pd.read_csv(args.output_dir / "synthesis" / "shared_exact_species.csv")
    require(len(shared) == expected["shared_exact_species"], "Shared-species count differs")
    print("All reference checkpoints passed.")


if __name__ == "__main__":
    main()
