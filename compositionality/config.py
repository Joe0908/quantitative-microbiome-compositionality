"""Fixed analysis configuration used in the completed project."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = PROJECT_ROOT / "data"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"

SEED = 531
N_SPLITS = 5
N_REPEATS = 10
PRIMARY_PREVALENCE = 0.05
LCPM_DETECTION_LIMIT = 1e-6

HGB_PARAMETERS = {
    "learning_rate": 0.05,
    "max_iter": 150,
    "max_leaf_nodes": 15,
    "min_samples_leaf": 10,
    "l2_regularization": 1.0,
    "early_stopping": False,
    "random_state": SEED,
}

LCPM_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41591-024-02963-2/MediaObjects/"
    "41591_2024_2963_MOESM3_ESM.xlsx"
)
LCPM_SHA256 = "2d9c0fa807fbcd85f97beee292d24551920d33bc76435c4ea578f5d90cc10282"

METACARDIS_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41591-022-01688-4/MediaObjects/"
    "41591_2022_1688_MOESM3_ESM.xlsx"
)
METACARDIS_SHA256 = "ce68279db4ce3c0c29a244ef6fb5ff712dcee3e6c3b5b33250174368b0c74248"

LCPM_FILENAME = "LCPM_supplementary_tables.xlsx"
METACARDIS_FILENAME = "MetaCardis_supplementary_tables.xlsx"

CORE_COVARIATES = [
    "age",
    "bmi",
    "bmi_missing",
    "male",
    "diabetes",
]

MEDICATION_COVARIATES = [
    "antidiabetic_drug",
    "antihypertensive_drug",
    "lipid_lowering_drug",
    "proton_pump_inhibitor",
]
