"""Part 02 — construct paired LCPM QMP and derived RMP matrices."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .common import BINOMIAL_DOT, close_rows, ensure_dir, require, write_json
from .config import DEFAULT_DATA_DIR, LCPM_FILENAME


DIAGNOSIS_LABELS = {
    "CTL": "No colonic lesions (control)",
    "ADE": "Adenoma / colorectal polyps",
    "CRC": "Colorectal cancer",
}


def run(data_dir: Path = DEFAULT_DATA_DIR) -> dict[str, Path]:
    raw_path = data_dir / "raw" / LCPM_FILENAME
    processed_dir = ensure_dir(data_dir / "processed" / "lcpm")

    metadata = pd.read_excel(raw_path, sheet_name="S1", header=0)
    metadata = metadata.rename(
        columns={
            "sample_ID": "participant_id",
            "Diagnosis": "diagnosis",
            "Colonoscopy": "colonoscopy_referral",
        }
    )
    metadata["participant_id"] = metadata["participant_id"].astype(str)
    metadata["diagnosis_label"] = metadata["diagnosis"].map(DIAGNOSIS_LABELS)
    metadata["progression_order"] = metadata["diagnosis"].map(
        {"CTL": 0, "ADE": 1, "CRC": 2}
    )

    qmp = pd.read_excel(raw_path, sheet_name="S14", header=0)
    qmp = qmp.rename(columns={qmp.columns[0]: "participant_id"})
    qmp["participant_id"] = qmp["participant_id"].astype(str)
    qmp = qmp.set_index("participant_id").apply(pd.to_numeric, errors="raise")
    qmp = qmp.loc[metadata["participant_id"]]
    rmp = close_rows(qmp)

    s6 = pd.read_excel(raw_path, sheet_name="S6", header=0)
    clean_names = s6["Species name"].dropna().astype(str).replace(
        {"Escherichia.Shigella.coli": "Escherichia-Shigella.coli"}
    )
    clean_set = set(clean_names)
    require(len(clean_set) == 336, "Expected 336 LCPM clean candidate labels")

    taxonomy = pd.DataFrame({"feature_id": qmp.columns})
    taxonomy["clean_candidate"] = taxonomy["feature_id"].isin(clean_set)
    taxonomy["strict_binomial"] = taxonomy["feature_id"].map(
        lambda value: bool(BINOMIAL_DOT.fullmatch(str(value)))
    )
    taxonomy["canonical_species"] = taxonomy["feature_id"].where(
        taxonomy["strict_binomial"]
    )
    taxonomy["canonical_species"] = taxonomy["canonical_species"].str.replace(
        ".", " ", regex=False
    )

    metadata = metadata.set_index("participant_id").loc[qmp.index].reset_index()
    metadata["total_load_cells_per_g"] = qmp.sum(axis=1).to_numpy()
    metadata["detected_features_n"] = (qmp > 0).sum(axis=1).to_numpy()
    metadata["derived_relative_sum"] = rmp.sum(axis=1).to_numpy()

    require(qmp.shape == (589, 676), f"Unexpected LCPM QMP shape: {qmp.shape}")
    require(metadata["diagnosis"].value_counts().to_dict() == {"ADE": 337, "CTL": 205, "CRC": 47}, "Unexpected LCPM diagnosis counts")
    require(np.allclose(rmp.sum(axis=1), 1.0), "LCPM RMP rows do not sum to one")

    paths = {
        "metadata": processed_dir / "metadata.csv.gz",
        "taxonomy": processed_dir / "taxonomy.csv.gz",
        "qmp": processed_dir / "qmp_cells_per_g.csv",
        "rmp": processed_dir / "rmp_row_closed.csv",
        "qc": processed_dir / "qc.json",
    }
    metadata.to_csv(paths["metadata"], index=False)
    taxonomy.to_csv(paths["taxonomy"], index=False)
    qmp.to_csv(paths["qmp"], index=True)
    rmp.to_csv(paths["rmp"], index=True)
    write_json(
        {
            "samples": 589,
            "features": 676,
            "clean_candidates": 336,
            "diagnosis_counts": metadata["diagnosis"].value_counts().to_dict(),
            "maximum_rmp_row_sum_error": float(np.abs(rmp.sum(axis=1) - 1).max()),
        },
        paths["qc"],
    )
    print(f"prepared LCPM: {qmp.shape[0]} samples × {qmp.shape[1]} features")
    return paths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    args = parser.parse_args()
    run(args.data_dir)


if __name__ == "__main__":
    main()
