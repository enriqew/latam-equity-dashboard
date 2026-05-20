"""Ingest PGx artefacts from the pgx-latam-atlas pipeline output.

Expected source: the portfolio repo's src/data/pgx/ directory (local copy)
or the published GitHub release artefacts from enriqew/pgx-latam-atlas.

Files fetched / copied:
  - drug_impact_summary.json
  - actionability_ranking.json

Usage:
    python -m ingest.pgx_artifacts [--source local|github]
"""
from __future__ import annotations

import argparse
import json
import logging
import shutil
from pathlib import Path

from ingest._common import REPO_ROOT, RAW_DIR, ensure_dirs, fetch_json, log, save_raw

# Local portfolio path (sibling repo or configured path)
PORTFOLIO_PGX_DIR = Path(r"C:\Users\ciker\PycharmProjects\data-dive-design-hub\src\data\pgx")

# GitHub raw artefact base URL (update tag as pipeline releases)
GITHUB_RAW_BASE = (
    "https://raw.githubusercontent.com/enriqew/pgx-latam-atlas/main/artifacts"
)

PGX_FILES = [
    "drug_impact_summary.json",
    "actionability_ranking.json",
]


def ingest_from_local() -> None:
    """Copy PGx JSON files from the local portfolio repo."""
    ensure_dirs()
    if not PORTFOLIO_PGX_DIR.exists():
        raise FileNotFoundError(
            f"Local portfolio PGx dir not found: {PORTFOLIO_PGX_DIR}\n"
            "Check the path or use --source github."
        )
    for fname in PGX_FILES:
        src = PORTFOLIO_PGX_DIR / fname
        if not src.exists():
            log.warning("File not found in local source: %s", src)
            continue
        dst = RAW_DIR / fname
        shutil.copy2(src, dst)
        log.info("Copied %s → %s", src, dst)


def ingest_from_github() -> None:
    """Fetch PGx JSON files from the GitHub raw URL."""
    ensure_dirs()
    for fname in PGX_FILES:
        url = f"{GITHUB_RAW_BASE}/{fname}"
        log.info("Fetching %s …", url)
        data = fetch_json(url)
        save_raw(data, fname)


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description="Ingest PGx artefacts")
    parser.add_argument(
        "--source",
        choices=["local", "github"],
        default="local",
        help="Where to fetch PGx files from (default: local portfolio repo)",
    )
    args = parser.parse_args(argv)

    if args.source == "local":
        ingest_from_local()
    else:
        ingest_from_github()

    log.info("PGx ingest complete. Files in %s", RAW_DIR)


if __name__ == "__main__":
    main()
