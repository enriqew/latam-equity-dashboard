"""Build per-country transplant time series for the LATAM equity dashboard.

Produces a list of country-level series objects sorted by country name,
each containing year-by-year transplant volume records.  European reference
countries (ESP, GBR) are included so the dashboard can render the access gap
visually over time.
"""
from __future__ import annotations

LATAM_COUNTRIES: dict[str, str] = {
    "MEX": "MXL",
    "PER": "PEL",
    "COL": "CLM",
}

REFERENCE_COUNTRIES: set[str] = {"ESP", "GBR"}

SERIES_FIELDS = (
    "year",
    "transplants_pmp",
    "total_transplants",
    "deceased_donors_pmp",
    "living_donors_pmp",
)


def build_timeseries(world_transplants: list[dict]) -> list[dict]:
    """
    Group transplant rows into per-country time series.

    Args:
        world_transplants: rows from world-transplants.json
    Returns:
        List of country objects with a nested ``series`` list, sorted by
        country_iso3.  Only rows with a non-null transplants_pmp are included.
    """
    target_iso3s = set(LATAM_COUNTRIES) | REFERENCE_COUNTRIES

    buckets: dict[str, dict] = {}
    for row in world_transplants:
        iso3 = row.get("country_iso3", "")
        if iso3 not in target_iso3s:
            continue
        if row.get("transplants_pmp") is None:
            continue

        if iso3 not in buckets:
            buckets[iso3] = {
                "country_iso3": iso3,
                "country_name": row.get("country_name", iso3),
                "population_code": LATAM_COUNTRIES.get(iso3),
                "is_reference": iso3 in REFERENCE_COUNTRIES,
                "series": [],
            }

        point = {field: row.get(field) for field in SERIES_FIELDS}
        buckets[iso3]["series"].append(point)

    for entry in buckets.values():
        entry["series"].sort(key=lambda p: p["year"] or 0)

    return sorted(buckets.values(), key=lambda e: e["country_iso3"])
