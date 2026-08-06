"""Part 01 — download and checksum the two public supplementary workbooks."""

from __future__ import annotations

import argparse
import shutil
import urllib.request
from pathlib import Path

from .common import ensure_dir, sha256sum
from .config import (
    DEFAULT_DATA_DIR,
    LCPM_FILENAME,
    LCPM_SHA256,
    LCPM_URL,
    METACARDIS_FILENAME,
    METACARDIS_SHA256,
    METACARDIS_URL,
)


SOURCES = {
    LCPM_FILENAME: (LCPM_URL, LCPM_SHA256),
    METACARDIS_FILENAME: (METACARDIS_URL, METACARDIS_SHA256),
}


def download_one(url: str, expected_sha256: str, destination: Path) -> Path:
    if destination.exists() and sha256sum(destination) == expected_sha256:
        print(f"verified existing: {destination}")
        return destination
    ensure_dir(destination.parent)
    temporary = destination.with_suffix(destination.suffix + ".partial")
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request) as response, temporary.open("wb") as output:
        shutil.copyfileobj(response, output)
    observed = sha256sum(temporary)
    if observed != expected_sha256:
        temporary.unlink(missing_ok=True)
        raise ValueError(
            f"Checksum mismatch for {url}: expected {expected_sha256}, observed {observed}"
        )
    temporary.replace(destination)
    print(f"downloaded and verified: {destination}")
    return destination


def run(data_dir: Path = DEFAULT_DATA_DIR) -> list[Path]:
    raw_dir = ensure_dir(data_dir / "raw")
    return [
        download_one(url, checksum, raw_dir / filename)
        for filename, (url, checksum) in SOURCES.items()
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    args = parser.parse_args()
    run(args.data_dir)


if __name__ == "__main__":
    main()

