"""Export pipeline: load processed data, compute equity index, write JSON artefacts.

Usage:
    python -m src.latam_equity.export

Output files (written to artifacts/):
    latam_equity_index.json          — 4 countries, 1000G+CPIC PGx (original)
    latam_equity_index_extended.json — 8 countries, gnomAD HWE PGx (consistent methodology)
    latam_transplant_timeseries.json — per-country transplant time series
    latam_equity_meta.json           — run metadata
"""
from __future__ import annotations

import dataclasses
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

from src.latam_equity.equity_index import build_equity_points
from src.latam_equity.equity_index_extended import build_equity_points_extended
from src.latam_equity.timeseries import build_timeseries

log = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
ARTIFACTS_DIR = REPO_ROOT / "artifacts"

TRANSPLANTS_PATH = DATA_DIR / "raw" / "world-transplants.json"
PGX_PATH = DATA_DIR / "raw" / "drug_impact_summary.json"
GNOMAD_PGX_PATH = DATA_DIR / "raw" / "gnomad_cyp3a5.json"

OUTPUT_EQUITY = ARTIFACTS_DIR / "latam_equity_index.json"
OUTPUT_EQUITY_EXTENDED = ARTIFACTS_DIR / "latam_equity_index_extended.json"
OUTPUT_TIMESERIES = ARTIFACTS_DIR / "latam_transplant_timeseries.json"
OUTPUT_META = ARTIFACTS_DIR / "latam_equity_meta.json"


def _load_json(path: Path) -> list[dict] | dict:
    if not path.exists():
        raise FileNotFoundError(
            f"Source file not found: {path}\n"
            "Run `make ingest` first to download raw data."
        )
    with path.open(encoding="utf-8") as fh:
        data = json.load(fh)
    return data if isinstance(data, (list, dict)) else []


def _load_list(path: Path) -> list[dict]:
    data = _load_json(path)
    return data if isinstance(data, list) else data.get("data", [])


def run() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)

    log.info("Loading source files …")
    transplants = _load_list(TRANSPLANTS_PATH)
    drug_impact = _load_list(PGX_PATH)
    gnomad_pgx = _load_json(GNOMAD_PGX_PATH)
    log.info("  transplants: %d rows", len(transplants))
    log.info("  drug_impact: %d rows", len(drug_impact))
    log.info("  gnomad_pgx: %d populations", len(gnomad_pgx.get("populations", [])) if isinstance(gnomad_pgx, dict) else 0)

    log.info("Computing equity index (4 countries, 1000G PGx) …")
    points = build_equity_points(transplants, drug_impact)

    records = [dataclasses.asdict(p) for p in points]
    OUTPUT_EQUITY.write_text(
        json.dumps(records, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    log.info("Wrote %s (%d records)", OUTPUT_EQUITY, len(records))

    log.info("Computing extended equity index (8 countries, gnomAD HWE PGx) …")
    points_ext = build_equity_points_extended(transplants, gnomad_pgx)
    records_ext = [dataclasses.asdict(p) for p in points_ext]
    OUTPUT_EQUITY_EXTENDED.write_text(
        json.dumps(records_ext, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    log.info("Wrote %s (%d records)", OUTPUT_EQUITY_EXTENDED, len(records_ext))

    log.info("Building transplant time series …")
    timeseries = build_timeseries(transplants)
    OUTPUT_TIMESERIES.write_text(
        json.dumps(timeseries, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    log.info(
        "Wrote %s (%d countries)",
        OUTPUT_TIMESERIES,
        len(timeseries),
    )

    meta = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_transplants": str(TRANSPLANTS_PATH.relative_to(REPO_ROOT)),
        "source_pgx": str(PGX_PATH.relative_to(REPO_ROOT)),
        "record_count": len(records),
        "populations": [r["population_code"] for r in records],
    }
    OUTPUT_META.write_text(
        json.dumps(meta, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    log.info("Wrote %s", OUTPUT_META)


if __name__ == "__main__":
    run()
