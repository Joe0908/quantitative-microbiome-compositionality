"""Part 03 — audit MetaCardis labels and construct QMP/row-closed matrices."""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import DEFAULT_DATA_DIR, METACARDIS_FILENAME


MEMBERSHIP_COLUMNS = {
    "HC275": "hc_member",
    "MMC372": "mmc_member",
    "UMCC222": "ummc_member",
    "IHD372": "ihd_member",
    "ACS112": "acs_member",
    "CIHD158": "cihd_member",
    "HF102": "hf_member",
}


def _assert_duplicate_rows_identical(frame: pd.DataFrame, name: str) -> None:
    value_columns = [column for column in frame.columns if column not in {"ID", "Status"}]
    for participant_id, positions in frame.groupby("ID", sort=False).indices.items():
        if len(positions) < 2:
            continue
        reference = frame.iloc[positions[0]][value_columns]
        for position in positions[1:]:
            if not reference.equals(frame.iloc[position][value_columns]):
                raise ValueError(
                    f"{name}: participant {participant_id} differs beyond Status"
                )


def _first_by_id(frame: pd.DataFrame) -> pd.DataFrame:
    first_rows = frame.loc[~frame["ID"].duplicated(keep="first")]
    return first_rows.drop(columns="Status").set_index("ID")


def _audit_stratum(row: pd.Series) -> str:
    if row["ihd_member"]:
        return "IHD"
    if row["hc_member"] and not row["mmc_member"] and not row["ummc_member"]:
        return "HC only"
    if row["mmc_member"]:
        return "MMC"
    if row["ummc_member"]:
        return "UMMC"
    return "Other"


def run(data_dir: Path = DEFAULT_DATA_DIR) -> dict[str, Path]:
    raw_path = data_dir / "raw" / METACARDIS_FILENAME

    global pd, np, BINOMIAL_SPACE, cag_key, close_rows, ensure_dir, require, write_json, yes_no_to_binary
    import numpy as np
    import pandas as pd
    from .common import (
        BINOMIAL_SPACE,
        cag_key,
        close_rows,
        ensure_dir,
        require,
        write_json,
        yes_no_to_binary,
    )
    from .xlsx_lite import LiteXlsx

    # Stream cached values directly from worksheet XML. The public workbook is
    # heavily formatted; this avoids loading an unnecessary full style model.
    with LiteXlsx(raw_path) as workbook:
        st9 = workbook.read_sheet("ST9", header=1)
        st10 = workbook.read_sheet("ST10", header=1)
        st14 = workbook.read_sheet("ST14", header=1)
        st5 = workbook.read_sheet("ST5", header=1, usecols_end=8)
    print("streamed required MetaCardis sheets", flush=True)
    processed_dir = ensure_dir(data_dir / "processed" / "metacardis")

    for frame, name in [(st9, "ST9"), (st10, "ST10"), (st14, "ST14")]:
        frame["ID"] = frame["ID"].astype(str)
        _assert_duplicate_rows_identical(frame, name)
        require(len(frame) == 1882, f"{name}: expected 1,882 published label rows")
        require(frame["ID"].nunique() == 1087, f"{name}: expected 1,087 unique IDs")
        print(f"audited {name}", flush=True)

    memberships = (
        st9.groupby("ID")["Status"]
        .agg(lambda values: tuple(dict.fromkeys(values.astype(str))))
        .rename("memberships")
        .to_frame()
    )
    memberships["published_label_rows"] = st9.groupby("ID").size()
    for status, column in MEMBERSHIP_COLUMNS.items():
        memberships[column] = memberships["memberships"].map(lambda values: status in values)
    memberships["published_analysis_memberships"] = memberships["memberships"].map("; ".join)
    memberships["non_overlapping_audit_stratum"] = memberships.apply(_audit_stratum, axis=1)

    metadata = _first_by_id(st9)
    drugs = _first_by_id(st14)
    taxonomy_meta = _first_by_id(st10.iloc[:, :5])
    metadata = memberships.join(metadata, how="left").join(
        taxonomy_meta.add_prefix("taxonomy__"), how="left"
    ).join(drugs.add_prefix("drug__"), how="left")

    feature_columns = st10.columns[5:].tolist()
    qmp_all = _first_by_id(st10[["ID", "Status", *feature_columns]])
    qmp_all = qmp_all.apply(pd.to_numeric, errors="coerce")
    complete = qmp_all.notna().all(axis=1) & (qmp_all.sum(axis=1) > 0)
    qmp = qmp_all.loc[complete]
    rmp = close_rows(qmp)

    metadata["quantitative_profile_available"] = metadata.index.isin(qmp.index)
    load = pd.to_numeric(metadata["taxonomy__Microbial load"], errors="coerce")
    metadata["microbial_load_available"] = load.notna() & (load > 0)
    metadata["quantitative_load_ready"] = (
        metadata["quantitative_profile_available"]
        & metadata["microbial_load_available"]
    )

    taxonomy_rows = st5.loc[
        st5["Feature ID"].astype(str).str.contains(r"_CAG", regex=True, na=False)
    ].copy()
    taxonomy_rows["cag_key"] = taxonomy_rows["Feature ID"].map(cag_key)
    matrix_dictionary = pd.DataFrame({"matrix_column": feature_columns})
    matrix_dictionary["cag_key"] = matrix_dictionary["matrix_column"].map(cag_key)
    taxonomy = matrix_dictionary.merge(taxonomy_rows, on="cag_key", how="left", validate="one_to_one")
    require(taxonomy["Feature ID"].notna().all(), "Failed to map every MetaCardis MGS feature to ST5 taxonomy")
    taxonomy["clean_candidate"] = (
        taxonomy["superkingdom"].eq("Bacteria")
        & ~taxonomy["species"].astype(str).str.lower().str.startswith("unclassified")
    )
    taxonomy["strict_binomial"] = taxonomy["species"].map(
        lambda value: bool(BINOMIAL_SPACE.fullmatch(str(value)))
    )
    # Conservative source-harmonization flag used by the completed project.
    # It removes generic ``sp.``/``bacterium`` labels before the stricter
    # two-token binomial and one-to-one rules are applied in Part 07.
    taxonomy["cross_cohort_exact_species_candidate"] = (
        taxonomy["clean_candidate"]
        & ~taxonomy["species"].astype(str).str.contains(
            r"\bsp\.|bacterium", case=False, regex=True
        )
    )

    rename = {
        "memberships": "published_membership_tuple",
        "Age (years)": "age",
        "BMI (kg/m²)": "bmi",
        "Gender": "gender",
        "Nationality": "nationality",
        "Diabetic status": "diabetic_status",
        "taxonomy__MGS count": "mgs_count",
        "taxonomy__Gene count": "gene_count",
        "taxonomy__Microbial load": "microbial_load_cells_per_g",
        "drug__Antidiabetic drugs": "antidiabetic_drug_raw",
        "drug__Anti-hypertensive drugs": "antihypertensive_drug_raw",
        "drug__Lipid lowering drugs": "lipid_lowering_drug_raw",
        "drug__Proton pump inhibitor": "proton_pump_inhibitor_raw",
    }
    metadata = metadata.rename(columns=rename)
    for column in ["age", "bmi", "microbial_load_cells_per_g"]:
        metadata[column] = pd.to_numeric(metadata[column], errors="coerce")
    metadata["male"] = metadata["gender"].astype(str).str.strip().str.lower().eq("male").astype(float)
    metadata["diabetes"] = yes_no_to_binary(metadata["diabetic_status"])
    for output, source in [
        ("antidiabetic_drug", "antidiabetic_drug_raw"),
        ("antihypertensive_drug", "antihypertensive_drug_raw"),
        ("lipid_lowering_drug", "lipid_lowering_drug_raw"),
        ("proton_pump_inhibitor", "proton_pump_inhibitor_raw"),
    ]:
        metadata[output] = yes_no_to_binary(metadata[source])

    require(qmp.shape == (994, 729), f"Unexpected MetaCardis QMP shape: {qmp.shape}")
    require(int(taxonomy["clean_candidate"].sum()) == 416, "Expected 416 clean bacterial MetaCardis candidates")
    require(
        np.allclose(rmp.sum(axis=1), 1.0),
        "MetaCardis QMP-derived row-closed rows do not sum to one",
    )
    primary_ids = metadata.index[
        metadata["quantitative_profile_available"]
        & (metadata["ihd_member"] | metadata["mmc_member"])
    ]
    primary = metadata.loc[primary_ids]
    require(not (primary["ihd_member"] & primary["mmc_member"]).any(), "IHD372 and MMC372 must be ID-disjoint")
    require(int(primary["ihd_member"].sum()) == 303, "Expected 303 IHD profiles")
    require(int(primary["mmc_member"].sum()) == 369, "Expected 369 MMC profiles")

    paths = {
        "metadata": processed_dir / "metadata.csv.gz",
        "taxonomy": processed_dir / "taxonomy.csv.gz",
        "qmp": processed_dir / "qmp_index.csv",
        "row_closed": processed_dir / "row_closed.csv",
        "qc": processed_dir / "qc.json",
    }
    metadata.reset_index(names="participant_id").to_csv(paths["metadata"], index=False)
    taxonomy.to_csv(paths["taxonomy"], index=False)
    qmp.to_csv(paths["qmp"], index=True)
    rmp.to_csv(paths["row_closed"], index=True)
    write_json(
        {
            "published_label_rows": 1882,
            "unique_participant_ids": 1087,
            "complete_quantitative_profiles": int(len(qmp)),
            "features": int(qmp.shape[1]),
            "clean_candidates": int(taxonomy["clean_candidate"].sum()),
            "primary_ihd_profiles": int(primary["ihd_member"].sum()),
            "primary_mmc_profiles": int(primary["mmc_member"].sum()),
            "maximum_row_closed_sum_error": float(
                np.abs(rmp.sum(axis=1) - 1).max()
            ),
        },
        paths["qc"],
    )
    print(f"prepared MetaCardis: {qmp.shape[0]} profiles × {qmp.shape[1]} MGS features")
    return paths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    args = parser.parse_args()
    run(args.data_dir)


if __name__ == "__main__":
    main()
