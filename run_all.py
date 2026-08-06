"""Run the complete quantitative-versus-compositional analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

from compositionality.config import DEFAULT_DATA_DIR, DEFAULT_OUTPUT_DIR
from compositionality.part01_download_data import run as download
from compositionality.part02_prepare_lcpm import run as prepare_lcpm
from compositionality.part03_prepare_metacardis import run as prepare_metacardis
from compositionality.part04_lcpm_associations import run as associate_lcpm
from compositionality.part05_metacardis_associations import run as associate_metacardis
from compositionality.part06_train_models import run as train_models
from compositionality.part07_synthesize import run as synthesize


def _part(number: int, label: str) -> None:
    print(f"\nPART {number:02d} — {label}\n{'=' * 72}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reproduce the LCPM + MetaCardis compositionality project."
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Use already-downloaded workbooks in DATA_DIR/raw.",
    )
    parser.add_argument(
        "--microbiome-only",
        action="store_true",
        help="Skip optional MetaCardis clinical and combined prediction models.",
    )
    args = parser.parse_args()

    if not args.skip_download:
        _part(1, "download and verify public inputs")
        download(args.data_dir)
    else:
        print("\nPART 01 — skipped by request", flush=True)

    _part(2, "prepare Galazzo/LCPM")
    prepare_lcpm(args.data_dir)
    _part(3, "prepare and audit MetaCardis")
    prepare_metacardis(args.data_dir)
    _part(4, "LCPM differential associations")
    associate_lcpm(args.data_dir, args.output_dir)
    _part(5, "MetaCardis hurdle and CLR associations")
    associate_metacardis(args.data_dir, args.output_dir)
    _part(6, "repeated cross-validation")
    train_models(
        args.data_dir,
        args.output_dir,
        include_metacardis_clinical=not args.microbiome_only,
    )
    _part(7, "guarded cross-cohort synthesis")
    synthesize(args.output_dir)
    print(f"\nComplete. Results are in: {args.output_dir.resolve()}", flush=True)


if __name__ == "__main__":
    main()
