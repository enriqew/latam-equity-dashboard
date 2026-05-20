"""Shared utilities for ingest modules: HTTP helpers, caching, schema validation."""
from __future__ import annotations

import hashlib
import json
import logging
import time
from pathlib import Path
from typing import Any

import requests

log = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = REPO_ROOT / "data" / "raw"
SCHEMAS_DIR = REPO_ROOT / "schemas"

# Default request settings
DEFAULT_TIMEOUT = 30
DEFAULT_RETRIES = 3
DEFAULT_BACKOFF = 2.0  # seconds, doubled each retry


def ensure_dirs() -> None:
    """Create data/raw if it doesn't exist."""
    RAW_DIR.mkdir(parents=True, exist_ok=True)


def fetch_json(
    url: str,
    *,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    backoff: float = DEFAULT_BACKOFF,
    headers: dict[str, str] | None = None,
) -> Any:
    """
    GET a URL and return the parsed JSON.  Retries with exponential backoff.

    Args:
        url: Full URL to fetch.
        timeout: Per-request timeout in seconds.
        retries: Number of retry attempts after initial failure.
        backoff: Base backoff time (doubled each attempt).
        headers: Extra HTTP headers.

    Returns:
        Parsed JSON (dict or list).

    Raises:
        requests.HTTPError: On non-2xx after all retries.
        requests.Timeout: If every attempt times out.
    """
    _headers = {"Accept": "application/json"}
    if headers:
        _headers.update(headers)

    last_exc: Exception | None = None
    for attempt in range(retries + 1):
        try:
            resp = requests.get(url, headers=_headers, timeout=timeout)
            resp.raise_for_status()
            return resp.json()
        except (requests.HTTPError, requests.Timeout, requests.ConnectionError) as exc:
            last_exc = exc
            if attempt < retries:
                wait = backoff * (2 ** attempt)
                log.warning("Attempt %d/%d failed (%s). Retrying in %.1fs …", attempt + 1, retries + 1, exc, wait)
                time.sleep(wait)

    raise RuntimeError(f"All {retries + 1} attempts failed for {url}") from last_exc


def save_raw(data: Any, filename: str) -> Path:
    """
    Write data as pretty-printed JSON to data/raw/<filename>.

    Returns the written path.
    """
    ensure_dirs()
    path = RAW_DIR / filename
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    log.info("Saved %s (%d bytes)", path, path.stat().st_size)
    return path


def load_raw(filename: str) -> Any:
    """Read and parse data/raw/<filename>."""
    path = RAW_DIR / filename
    if not path.exists():
        raise FileNotFoundError(f"Raw file not found: {path}. Run ingest first.")
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def checksum(path: Path) -> str:
    """Return SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()
