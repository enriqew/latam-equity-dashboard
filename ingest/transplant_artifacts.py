"""Ingest world transplant data from the transplant-atlas pipeline output.

Expected source: the portfolio repo's src/data/transplant-atlas/ directory (local copy)
or the published GitHub release artefacts from enriqew/transplant-atlas.

Files fetched / copied:
  - world-transplants.json

Usage:
    python -m ingest.transplant_artifacts [--source local|github]
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil
from pathlib import Path

from ingest._common import REPO_ROOT, RAW_DIR, ensure_dirs, fetch_json, log, save_raw

# Local portfolio path. Defaults to a sibling checkout of the portfolio repo;
# override with PORTFOLIO_DIR when it lives somewhere else.
PORTFOLIO_DIR = Path(os.environ.get("PORTFOLIO_DIR", REPO_ROOT.parent / "data-dive-design-hub"))
PORTFOLIO_TRANSPLANT_DIR = PORTFOLIO_DIR / "src" / "data" / "transplant-atlas"

# GitHub raw artefact base URL
GITHUB_RAW_BASE = (
    "https://raw.githubusercontent.com/enriqew/transplant-atlas/main/data/exports"
)

TRANSPLANT_FILES = [
    "world-transplants.json",
]


def ingest_from_local() -> None:
    """Copy transplant JSON files from the local portfolio repo."""
    ensure_dirs()
    if not PORTFOLIO_TRANSPLANT_DIR.exists():
        raise FileNotFoundError(
            f"Local portfolio transplant dir not found: {PORTFOLIO_TRANSPLANT_DIR}\n"
            "Check the path or use --source github."
        )
    for fname in TRANSPLANT_FILES:
        src = PORTFOLIO_TRANSPLANT_DIR / fname
        if not src.exists():
            log.warning("File not found in local source: %s", src)
            continue
        dst = RAW_DIR / fname
        shutil.copy2(src, dst)
        log.info("Copied %s → %s", src, dst)


def ingest_from_github() -> None:
    """Fetch transplant JSON files from the GitHub raw URL."""
    ensure_dirs()
    for fname in TRANSPLANT_FILES:
        url = f"{GITHUB_RAW_BASE}/{fname}"
        log.info("Fetching %s …", url)
        data = fetch_json(url)
        save_raw(data, fname)


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description="Ingest transplant artefacts")
    parser.add_argument(
        "--source",
        choices=["local", "github"],
        default="local",
        help="Where to fetch transplant files from (default: local portfolio repo)",
    )
    args = parser.parse_args(argv)

    if args.source == "local":
        ingest_from_local()
    else:
        ingest_from_github()

    log.info("Transplant ingest complete. Files in %s", RAW_DIR)


if __name__ == "__main__":
    main()
