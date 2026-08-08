"""Part 04 — LCPM QMP, row-closed, and CLR association analyses."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .common import (
    bh_fdr,
    close_rows,
    ensure_dir,
    fixed_clr,
    multiplicative_clr,
    pooled_prevalence_filter,
    prevalence_filter,
    read_matrix,
    require,
    write_json,
)
from .config import (
    DEFAULT_DATA_DIR,
    DEFAULT_OUTPUT_DIR,
    LCPM_DETECTION_LIMIT,
    PRIMARY_PREVALENCE,
)


def _pairwise_results(
    component: str,
    representation: str,
    matrix: pd.DataFrame,
    labels: pd.Series,
    filter_specification: str,
    primary_analysis: bool,
    clr_zero_replacement: str = "not_applicable",
) -> pd.DataFrame:
    ctl_index = labels.index[labels == "CTL"]
    crc_index = labels.index[labels == "CRC"]
    rows: list[dict[str, object]] = []
    for feature in matrix.columns:
        ctl = matrix.loc[ctl_index, feature].to_numpy(dtype=float)
        crc = matrix.loc[crc_index, feature].to_numpy(dtype=float)
        test = stats.mannwhitneyu(ctl, crc, alternative="two-sided", method="asymptotic")
        effect = 1.0 - (2.0 * float(test.statistic) / (len(ctl) * len(crc)))
        rows.append(
            {
                "cohort": "LCPM",
                "contrast": "CRC vs CTL",
                "filter_specification": filter_specification,
                "primary_analysis": primary_analysis,
                "component": component,
                "representation": representation,
                "clr_zero_replacement": clr_zero_replacement,
                "feature": feature,
                "n_ctl": len(ctl),
                "n_crc": len(crc),
                "u_ctl": float(test.statistic),
                "p_value": float(test.pvalue),
                "effect": effect,
                "effect_definition": "rank-biserial CRC minus CTL",
                "median_ctl": float(np.median(ctl)),
                "median_crc": float(np.median(crc)),
                "median_difference_crc_minus_ctl": float(np.median(crc) - np.median(ctl)),
            }
        )
    result = pd.DataFrame(rows)
    result["q_value_bh"] = bh_fdr(result["p_value"])
    result["significant_q_lt_0_05"] = result["q_value_bh"] < 0.05
    return result.sort_values(["q_value_bh", "p_value", "feature"])


def _global_results(
    qmp: pd.DataFrame,
    rmp: pd.DataFrame,
    labels: pd.Series,
    candidates: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    groups = ["CTL", "ADE", "CRC"]
    for threshold in (0.01, 0.05, 0.10):
        features = prevalence_filter(
            rmp[candidates], labels, threshold, LCPM_DETECTION_LIMIT
        )
        clr = fixed_clr(rmp[features])
        specifications = [
            ("QMP", "Kruskal-Wallis", qmp[features]),
            ("Row-closed", "Kruskal-Wallis", rmp[features]),
            ("CLR", "one-way ANOVA", clr),
            ("CLR", "Kruskal-Wallis", clr),
        ]
        for representation, test_name, matrix in specifications:
            local_rows: list[dict[str, object]] = []
            for feature in features:
                arrays = [
                    matrix.loc[labels.index[labels == group], feature].to_numpy(dtype=float)
                    for group in groups
                ]
                if test_name == "one-way ANOVA":
                    test = stats.f_oneway(*arrays)
                else:
                    test = stats.kruskal(*arrays)
                local_rows.append(
                    {
                        "cohort": "LCPM",
                        "contrast": "CTL vs ADE vs CRC",
                        "filter_specification": "source_aligned_group_union",
                        "prevalence_threshold": threshold,
                        "representation": representation,
                        "test": test_name,
                        "feature": feature,
                        "statistic": float(test.statistic),
                        "p_value": float(test.pvalue),
                        "tested_features": len(features),
                    }
                )
            local = pd.DataFrame(local_rows)
            local["q_value_bh"] = bh_fdr(local["p_value"])
            local["significant_q_lt_0_05"] = local["q_value_bh"] < 0.05
            rows.extend(local.to_dict("records"))
    return pd.DataFrame(rows)


def run(
    data_dir: Path = DEFAULT_DATA_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> dict[str, Path]:
    processed = data_dir / "processed" / "lcpm"
    results_dir = ensure_dir(output_dir / "associations")
    metadata = pd.read_csv(processed / "metadata.csv.gz").set_index("participant_id")
    taxonomy = pd.read_csv(processed / "taxonomy.csv.gz")
    qmp = read_matrix(processed / "qmp_cells_per_g.csv")
    # Re-close the QMP matrix in memory. This preserves exact structural ties
    # that can be perturbed by text serialization at machine precision.
    rmp = close_rows(qmp)
    candidates = taxonomy.loc[taxonomy["clean_candidate"], "feature_id"].tolist()

    primary_ids = metadata.index[metadata["diagnosis"].isin(["CTL", "CRC"])]
    labels = metadata.loc[primary_ids, "diagnosis"]
    qmp_primary = qmp.loc[primary_ids, candidates]
    rmp_primary = rmp.loc[primary_ids, candidates]
    # The source-aligned analysis used the union of taxa detected in at least
    # 5% of any diagnosis group.  Because that rule uses outcome labels, it is
    # retained only as a sensitivity analysis.
    source_features = prevalence_filter(
        rmp[candidates],
        metadata["diagnosis"],
        PRIMARY_PREVALENCE,
        LCPM_DETECTION_LIMIT,
    )
    require(
        len(source_features) == 112,
        f"Expected 112 source-aligned LCPM features, observed {len(source_features)}",
    )

    # Primary inference uses a pooled prevalence calculation over the 252
    # participants in the contrast.  No diagnosis labels enter this filter.
    pooled_features = pooled_prevalence_filter(
        rmp_primary[candidates],
        PRIMARY_PREVALENCE,
        LCPM_DETECTION_LIMIT,
    )
    require(
        len(pooled_features) == 93,
        f"Expected 93 pooled outcome-blind LCPM features, observed {len(pooled_features)}",
    )

    pooled_clr_minimum = fixed_clr(rmp_primary[pooled_features])
    pooled_clr_multiplicative = multiplicative_clr(rmp_primary[pooled_features])
    source_clr_minimum = fixed_clr(rmp[source_features]).loc[primary_ids]
    source_clr_multiplicative = multiplicative_clr(rmp_primary[source_features])

    pairwise = pd.concat(
        [
            _pairwise_results(
                "qmp", "QMP", qmp_primary[pooled_features], labels,
                "pooled_outcome_blind", True,
            ),
            _pairwise_results(
                "row_closed", "Row-closed", rmp_primary[pooled_features], labels,
                "pooled_outcome_blind", True,
            ),
            _pairwise_results(
                "clr_minimum_positive", "CLR", pooled_clr_minimum, labels,
                "pooled_outcome_blind", True, "minimum_positive",
            ),
            _pairwise_results(
                "clr_multiplicative", "CLR", pooled_clr_multiplicative, labels,
                "pooled_outcome_blind", True, "multiplicative_delta_1_over_p_squared",
            ),
            _pairwise_results(
                "qmp", "QMP", qmp_primary[source_features], labels,
                "source_aligned_group_union", False,
            ),
            _pairwise_results(
                "row_closed", "Row-closed", rmp_primary[source_features], labels,
                "source_aligned_group_union", False,
            ),
            _pairwise_results(
                "clr_minimum_positive", "CLR", source_clr_minimum, labels,
                "source_aligned_group_union", False, "minimum_positive",
            ),
            _pairwise_results(
                "clr_multiplicative", "CLR", source_clr_multiplicative, labels,
                "source_aligned_group_union", False,
                "multiplicative_delta_1_over_p_squared",
            ),
        ],
        ignore_index=True,
    )
    global_da = _global_results(qmp, rmp, metadata["diagnosis"], candidates)

    pairwise_path = results_dir / "lcpm_crc_vs_ctl_associations.csv"
    global_path = results_dir / "lcpm_global_prevalence_sensitivity.csv"
    summary_path = results_dir / "lcpm_association_summary.json"
    pairwise.to_csv(pairwise_path, index=False)
    global_da.to_csv(global_path, index=False)
    counts = {}
    for specification, group in pairwise.groupby("filter_specification", sort=False):
        counts[specification] = (
            group.groupby("component")["significant_q_lt_0_05"]
            .sum()
            .astype(int)
            .to_dict()
        )
    write_json(
        {
            "contrast": "CRC vs CTL",
            "samples": int(len(primary_ids)),
            "crc": int((labels == "CRC").sum()),
            "ctl": int((labels == "CTL").sum()),
            "primary_filter": "pooled_outcome_blind",
            "primary_tested_features": len(pooled_features),
            "source_aligned_tested_features": len(source_features),
            "significant_calls": counts,
            "covariate_adjustment": "none; participant-level confounders were unavailable",
        },
        summary_path,
    )
    print(
        "LCPM associations: "
        f"{len(pooled_features)} pooled outcome-blind taxa and "
        f"{len(source_features)} source-aligned taxa; q<0.05 calls {counts}"
    )
    return {"pairwise": pairwise_path, "global": global_path, "summary": summary_path}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    run(args.data_dir, args.output_dir)


if __name__ == "__main__":
    main()
